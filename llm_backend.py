"""
llm_backend.py  —  Unified LLM backend for CoDHy
==================================================
Supports three runtime modes:
 
  • LOCAL      – model loaded from disk via HuggingFace transformers (GPU compute node,
                 no internet required at inference time).
                 Model must be pre-downloaded on the login node first.
 
  • HF_API     – calls the model via HuggingFace Inference Providers
                 (huggingface_hub.InferenceClient), routed to whichever
                 third-party provider (Together, Fireworks, SambaNova, etc.)
                 or HF's own "hf-inference" serves that model. No local GPU
                 or disk space needed. Requires HF_TOKEN.
                 NOTE: gating status and provider availability should be
                 verified on each model's HF page — MODEL_REGISTRY's
                 `gated` flag has been wrong before for at least one entry;
                 treat it as a starting guess, not a guarantee.
 
  • OLLAMA     – model served by a local Ollama daemon (legacy / dev mode).

Model registry
--------------------------
Set LLM_MODEL_KEY in secrets.env to one of the keys below.
The matching entry defines the repo, a (best-effort, verify yourself) gated
flag, VRAM estimate (only relevant for LOCAL), the HF Inference Provider to
pin (for HF_API — "auto" if unconfirmed), and notes.

KEY                     HF repo-id                                            Notes
----------------------  ----------------------------------------------------  -------
jsl_medllama_3_8b       johnsnowlabs/JSL-MedLlama-3-8B-v2.0                   biomedical; confirmed has an Inference Provider
bio_llama_3_8b          ContactDoctor/Bio-Medical-Llama-3-8B-CoT-012025       biomedical; gated (confirmed); no free Inference Provider as of last check
medgemma_27b_text_it    google/medgemma-27b-text-it                          biomedical; confirmed has an Inference Provider; needs transformers>=4.50.0 if run locally
meditron_7b             epfl-llm/meditron-7b                                  biomedical; NOT instruction-tuned (base model, no chat template) — poor fit for zero-shot JSON prompting
biom_mistral_7b         BioMistral/BioMistral-7B                              biomedical
pmc_llama_13b           axiong/PMC_LLaMA_13B                                  biomedical; gated
llama3_8b               meta-llama/Llama-3.1-8B-Instruct                      generic; gated (confirmed) — unclear if gating also blocks Inference Providers access, verify
deepseek_r1             deepseek-ai/DeepSeek-R1                               generic/reasoning; confirmed example model in HF's own Inference Providers docs; uses <think> reasoning blocks by design
qwen3_235b_instruct     Qwen/Qwen3-235B-A22B-Instruct-2507                    generic; MoE, 235B total/22B active params — must use hf_api, far too large for local

Usage in secrets.env
---------------------
  LLM_BACKEND=hf_api        # "local", "hf_api", or "ollama" (all three implemented)
  LLM_MODEL_KEY=medgemma_27b_text_it
  HF_TOKEN=hf_...            # required for hf_api; also needed for gated models under local
  HF_PROVIDER=auto           # or a specific provider name, e.g. "together", "sambanova" —
                             # pin this once you've checked which provider(s) actually list
                             # your chosen model, for run-to-run reproducibility
  LLM_LOCAL_DIR=/path/to/models   # only used when LLM_BACKEND=local
  OLLAMA_MODEL=llama3.1      # only used when LLM_BACKEND=ollama
  OLLAMA_BASE_URL=http://localhost:11434/api/generate
"""

import os
import json
import re
import logging
from typing import Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Model registry
# ---------------------------------------------------------------------------
# Switched from positional tuples to dicts -- the tuple had been extended
# three times already (gated -> +vram -> +notes -> +provider) and kept
# risking a silent positional mismatch. Dicts are self-documenting at each
# call site instead.
#
# call_style distinguishes two DIFFERENT, both real, both user-verified
# InferenceClient usage patterns (confirmed directly from each model's own
# HF page snippet, actually run successfully by the user):
#   "chat"           -> InferenceClient(api_key=...)  [no provider kwarg]
#                       client.chat.completions.create(
#                           model=f"{repo_id}:{provider}", messages=[...])
#   "text_generation"-> InferenceClient(provider=..., api_key=...)
#                       client.text_generation(prompt_string, model=repo_id,
#                                               max_new_tokens=...)
# These are NOT interchangeable -- do not assume one works for a model
# verified only under the other style.
MODEL_REGISTRY = {
    "qwen3_235b_instruct": {
        "repo_id": "Qwen/Qwen3-235B-A22B-Instruct-2507",
        "gated": None,  # not verified
        "provider": "novita",
        "call_style": "chat",
        "max_tokens_tested": 32768,  # user confirmed this value WORKED; true ceiling may be higher, untested
        "notes": (
            "Generic. Provider=novita, confirmed working by user at max_tokens=32768. "
            "'-Instruct' (not '-Thinking') name suggests no <think> block, but this is "
            "an unverified naming inference, not a confirmed fact."
        ),
    },
    "llama3_8b": {
        "repo_id": "meta-llama/Llama-3.1-8B-Instruct",
        "gated": True,  # gated on HF directly; confirmed separately that Inference Providers access still worked for user
        "provider": "novita",
        "call_style": "chat",
        "max_tokens_tested": 32768,  # user confirmed this value worked
        "notes": "Generic. Gated on HF, but user confirmed provider=novita access works with HF_TOKEN. max_tokens=32768 tested OK.",
    },
    "biomistral_7b_slerp": {
        "repo_id": "BioMistral/BioMistral-7B-SLERP",
        "gated": None,  # not verified
        "provider": "featherless-ai",
        "call_style": "chat",
        "max_tokens_tested": 4096,  # user-stated HARD ceiling for this model, not just "tested value" -- do not exceed
        "notes": (
            "Biomedical. Provider=featherless-ai. User states 4096 is this model's "
            "MAXIMUM token size (a real constraint, not just an untested ceiling) — "
            "requests above this will likely fail or be truncated by the provider."
        ),
    },
    "jsl_medllama_3_8b": {
        "repo_id": "johnsnowlabs/JSL-MedLlama-3-8B-v2.0",
        "gated": False,
        "provider": "featherless-ai",
        "call_style": "text_generation",  # NOTE: different call pattern, see module docstring
        "max_tokens_tested": 8000,  # user says "something around 8000" -- approximate, not exact
        "notes": (
            "Biomedical. Provider=featherless-ai. Uses client.text_generation(), NOT "
            "chat.completions.create() -- confirmed from the model's own HF page snippet. "
            "User says max token size is 'around 8000' -- approximate figure, not exact; "
            "verify the precise limit yourself if you need to rely on it exactly."
        ),
    },
}


def list_models() -> None:
    """Print the model registry in a readable format."""
    print(f"\n{'KEY':<22} {'REPO-ID':<45} {'GATED':<11} {'PROVIDER':<15} {'CALL STYLE':<17} {'MAX TOK':<8} NOTES")
    print("-" * 160)
    for key, info in MODEL_REGISTRY.items():
        gated = info["gated"]
        g = "UNVERIFIED" if gated is None else ("yes" if gated else "no")
        print(f"{key:<22} {info['repo_id']:<45} {g:<11} {info['provider']:<15} "
              f"{info['call_style']:<17} {info['max_tokens_tested']:<8} {info['notes']}")
    print()


# ---------------------------------------------------------------------------
# Backend implementations
# ---------------------------------------------------------------------------

class _OllamaBackend:
    """Thin wrapper around a local Ollama daemon (unchanged from original code)."""

    def __init__(self, model: str, base_url: str):
        import requests as _requests
        self._requests = _requests
        self.model = model
        self.api_url = base_url

    def generate(self, system_prompt: str, user_prompt: str,
                 temperature: float = 0.0, max_tokens: int = 3072,
                 json_mode: bool = False) -> str:
        payload = {
            "model": self.model,
            "prompt": system_prompt + "\n\n" + user_prompt,
            "stream": False,
            "options": {"temperature": temperature, "num_predict": max_tokens},
        }
        if json_mode:
            payload["format"] = "json"

        for attempt in range(3):
            try:
                r = self._requests.post(self.api_url, json=payload, timeout=300)
                if r.status_code == 200:
                    return r.json().get("response", "")
                logger.warning("Ollama HTTP %s on attempt %d", r.status_code, attempt + 1)
            except Exception as exc:
                logger.warning("Ollama connection error (attempt %d): %s", attempt + 1, exc)
        return ""

    def check_connection(self) -> bool:
        try:
            self._requests.get(self.api_url.replace("/generate", "/tags"), timeout=5)
            return True
        except Exception:
            return False


class _HFInferenceBackend:
    """
    Calls a model via HuggingFace's Inference Providers system, using
    huggingface_hub.InferenceClient. No local GPU or disk space needed.

    Supports TWO DIFFERENT, both real, both user-verified calling
    conventions, selected via `call_style` (do not assume these are
    interchangeable — each was confirmed by actually running the exact
    snippet from the model's own HF page):

      call_style="chat":
        InferenceClient(api_key=...)                     [NO provider kwarg here]
        client.chat.completions.create(
            model=f"{repo_id}:{provider}",                [provider embedded in model string]
            messages=[{"role": "user", "content": "..."}],  [SINGLE user message only --
                                                              system_prompt+user_prompt concatenated]
            max_tokens=..., temperature=...,
        )
        Confirmed working by user for: Qwen3-235B-A22B-Instruct-2507 (provider=novita),
        Llama-3.1-8B-Instruct (provider=novita), BioMistral-7B-SLERP (provider=featherless-ai).
        NOTE: an earlier version of this code used a separate "system" role message.
        None of the 3 verified snippets actually tested that -- all 3 used a single
        "user" message. A real error from featherless-ai ("Conversation roles must
        alternate user/assistant/...") confirmed the separate-system-role assumption
        was wrong, at least for that provider. Reverted to the only pattern with
        actual verified evidence: one "user" message, system+user text concatenated.

      call_style="text_generation":
        InferenceClient(provider=..., api_key=...)        [provider kwarg HERE instead]
        client.text_generation(
            prompt_string,                                [single string, NOT messages]
            model=repo_id,                                 [NO ":provider" suffix here]
            max_new_tokens=...,
        )
        Confirmed working by user for: JSL-MedLlama-3-8B-v2.0 (provider=featherless-ai).
        NOTE: this call style has no separate system/user role — system_prompt and
        user_prompt are concatenated into one string below. This is a plain
        completion, NOT officially chat-templated for this model — a design
        choice made here, not a verified best practice. If output quality
        seems off for this specific model, that concatenation is the first
        thing to reconsider.

    max_tokens_cap, if set, clamps requested max_tokens down to a per-model
    ceiling (BioMistral-7B-SLERP's 4096 is a real hard limit per the user;
    others are "tested OK at this value", not necessarily the true ceiling —
    see MODEL_REGISTRY notes per model for which is which).

    NOT verified (flagging rather than asserting):
      - Whether `response_format` (OpenAI-style JSON mode) is honored by
        chat_style providers/models here — only attempted for call_style="chat",
        never attempted for call_style="text_generation" (text_generation has
        no such parameter in the verified snippet, so it isn't tried).
      - Whether gating on meta-llama/Llama-3.1-8B-Instruct would block a
        DIFFERENT HF_TOKEN than the one the user already tested this with.
    """

    def __init__(self, repo_id: str, token: str = "", provider: str = "auto",
                 call_style: str = "chat", max_tokens_cap: Optional[int] = None):
        self.repo_id = repo_id
        self.token = token or None
        self.provider = provider
        self.call_style = call_style
        self.max_tokens_cap = max_tokens_cap
        self._client = None

    def _get_client(self):
        if self._client is None:
            from huggingface_hub import InferenceClient
            if self.call_style == "text_generation":
                # Provider goes at construction time for this call style
                # (verified from the JSL-MedLlama snippet).
                self._client = InferenceClient(provider=self.provider, api_key=self.token)
            else:
                # Provider is embedded in the model string at call time instead
                # (verified from the Qwen/Llama/BioMistral snippets) — no
                # provider kwarg here.
                self._client = InferenceClient(api_key=self.token)
        return self._client

    def _capped_max_tokens(self, requested: int) -> int:
        if self.max_tokens_cap is not None and requested > self.max_tokens_cap:
            logger.warning(
                "Requested max_tokens=%d exceeds this model's known cap of %d; "
                "clamping down.", requested, self.max_tokens_cap
            )
            return self.max_tokens_cap
        return requested

    def generate(self, system_prompt: str, user_prompt: str,
                 temperature: float = 0.0, max_tokens: int = 3072,
                 json_mode: bool = False) -> str:
        client = self._get_client()
        max_tokens = self._capped_max_tokens(max_tokens)

        last_exc = None
        for attempt in range(3):
            try:
                if self.call_style == "text_generation":
                    # No system/user role separation in this call style —
                    # concatenated. See class docstring caveat above.
                    prompt = f"{system_prompt}\n\n{user_prompt}"
                    result = client.text_generation(
                        prompt,
                        model=self.repo_id,
                        max_new_tokens=max_tokens,
                    )
                    return result or ""

                # call_style == "chat"
                # Fixed: none of the 3 verified snippets (Qwen, Llama-3.1,
                # BioMistral) ever used a separate "system" role — all 3
                # used a single "user" message. Adding a system role was an
                # unverified assumption on my part; a real 400 error
                # ("Conversation roles must alternate user/assistant/...")
                # confirmed it breaks at least featherless-ai's validator.
                # Since neither Qwen nor Llama's snippets tested a system
                # role either way, concatenating into one user message is
                # the only pattern with actual verified evidence behind it.
                messages = [
                    {"role": "user", "content": f"{system_prompt}\n\n{user_prompt}"},
                ]
                model_str = f"{self.repo_id}:{self.provider}"
                kwargs = {
                    "model": model_str,
                    "messages": messages,
                    "max_tokens": max_tokens,
                    "temperature": temperature,
                }
                if json_mode:
                    # Not verified to be honored by these providers/models —
                    # try it, fall back to plain generation if rejected.
                    try:
                        completion = client.chat.completions.create(
                            response_format={"type": "json_object"}, **kwargs
                        )
                        return completion.choices[0].message.content or ""
                    except TypeError:
                        pass  # kwarg itself rejected
                    except Exception as e:
                        logger.warning(
                            "response_format=json_object rejected (%s); retrying without it.", e
                        )

                completion = client.chat.completions.create(**kwargs)
                return completion.choices[0].message.content or ""

            except Exception as exc:
                last_exc = exc
                logger.warning("HF Inference error (attempt %d): %s", attempt + 1, exc)

        # Loud failure instead of silently returning "" -- an empty string
        # here previously surfaced several steps downstream as a confusing
        # "failed to parse JSON" message with no indication the real
        # problem was e.g. a 403 gating error or a provider outage.
        raise RuntimeError(
            f"HF Inference call failed after 3 attempts for model "
            f"'{self.repo_id}' (provider='{self.provider}', call_style='{self.call_style}'). "
            f"Last error: {last_exc}"
        )

    def check_connection(self) -> bool:
        try:
            client = self._get_client()
            if self.call_style == "text_generation":
                client.text_generation("ping", model=self.repo_id, max_new_tokens=1)
            else:
                client.chat.completions.create(
                    model=f"{self.repo_id}:{self.provider}",
                    messages=[{"role": "user", "content": "ping"}],
                    max_tokens=1,
                )
            return True
        except Exception as exc:
            logger.warning("HF Inference connection check failed: %s", exc)
            return False


class _LocalHFBackend:
    """
    Loads a model from a local directory using HuggingFace transformers.
    Lazy-loads on first call so import time stays fast.
    Works entirely offline — no internet required at inference time.
    """

    def __init__(self, model_dir: str, device: Optional[str] = None):
        self.model_dir = model_dir
        self._model = None
        self._tokenizer = None
        # Auto-detect device: prefer CUDA, fall back to CPU
        if device is None:
            try:
                import torch
                device = "cuda" if torch.cuda.is_available() else "cpu"
            except ImportError:
                device = "cpu"
        self.device = device
        logger.info("LocalHFBackend configured: dir=%s device=%s", model_dir, device)

    def _load(self):
        """Lazy-load the model + tokenizer once."""
        if self._model is not None:
            return
        try:
            import torch
            from transformers import AutoTokenizer, AutoModelForCausalLM
        except ImportError as e:
            raise RuntimeError(
                "transformers / torch not installed. "
                "Run: pip install transformers torch accelerate bitsandbytes"
            ) from e

        print(f"Loading model from {self.model_dir} onto {self.device} …")

        # 4-bit quantisation when on GPU — saves ~50% VRAM with minimal quality loss
        load_kwargs = {}
        if self.device.startswith("cuda"):
            try:
                from transformers import BitsAndBytesConfig
                load_kwargs["quantization_config"] = BitsAndBytesConfig(
                    load_in_4bit=True,
                    bnb_4bit_compute_dtype=torch.float16,
                )
                load_kwargs["device_map"] = "auto"
            except ImportError:
                load_kwargs["torch_dtype"] = torch.float16
                load_kwargs["device_map"] = "auto"

        self._tokenizer = AutoTokenizer.from_pretrained(self.model_dir, local_files_only=True)
        self._model = AutoModelForCausalLM.from_pretrained(
            self.model_dir, local_files_only=True, **load_kwargs
        )
        print("Model loaded successfully.")

    def generate(self, system_prompt: str, user_prompt: str,
                 temperature: float = 0.0, max_tokens: int = 3072,
                 json_mode: bool = False) -> str:
        """
        Follows the usage pattern verified directly from bio_llama_3_8b's
        HuggingFace model card: apply_chat_template(..., tokenize=True,
        return_dict=True, return_tensors="pt") -> model.generate(**inputs)
        -> tokenizer.decode(outputs[0][input_len:]).
        This replaces an earlier version that went through pipeline() and
        stripped the prompt via string-prefix matching, which is more
        fragile (depends on exact character-for-character round-tripping of
        the re-tokenized prompt string) than slicing by exact token count.

        NOTE on json_mode: vanilla transformers generate() has no built-in
        JSON-mode equivalent to Ollama's format="json" or OpenAI's
        response_format. json_mode is accepted here for interface
        consistency with the other backends but does not change behavior —
        true constrained decoding would need a grammar-enforcement library
        (e.g. outlines), which is a separate, bigger change.
        """
        self._load()
        import torch

        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ]

        # Robustness note (same reasoning as before): apply_chat_template
        # can raise for a model whose tokenizer defines no chat_template.
        # Verified working for bio_llama_3_8b directly from its model card,
        # but kept defensive in case LLM_MODEL_KEY changes to a model
        # without one — falls back rather than crashing.
        try:
            inputs = self._tokenizer.apply_chat_template(
                messages,
                add_generation_prompt=True,
                tokenize=True,
                return_dict=True,
                return_tensors="pt",
            ).to(self._model.device)
        except Exception as e:
            print(
                f"⚠️  WARNING: this model's tokenizer has no usable chat_template "
                f"({e}). Falling back to a plain string prompt, which may not "
                f"match what this model was actually trained on — check output "
                f"quality carefully."
            )
            full_prompt = (
                f"<|system|>\n{system_prompt}\n"
                f"<|user|>\n{user_prompt}\n"
                f"<|assistant|>\n"
            )
            inputs = self._tokenizer(full_prompt, return_tensors="pt").to(self._model.device)

        gen_kwargs = {
            "max_new_tokens": max_tokens,
            "do_sample": temperature > 0,
        }
        if temperature > 0:
            gen_kwargs["temperature"] = temperature
        gen_kwargs["pad_token_id"] = (
            self._tokenizer.pad_token_id
            if self._tokenizer.pad_token_id is not None
            else self._tokenizer.eos_token_id
        )

        with torch.no_grad():
            outputs = self._model.generate(**inputs, **gen_kwargs)

        # Slice off the prompt by exact token count (not string matching),
        # then decode only the newly generated tokens.
        input_len = inputs["input_ids"].shape[-1]
        new_tokens = outputs[0][input_len:]
        generated = self._tokenizer.decode(new_tokens, skip_special_tokens=True)
        return generated.strip()

    def check_connection(self) -> bool:
        return os.path.isdir(self.model_dir)


# ---------------------------------------------------------------------------
# Public factory
# ---------------------------------------------------------------------------

class LLMBackend:
    """
    Public interface used by DrugScreeningAgent and ValidationAgent.

    Instantiate once, pass around:
        backend = LLMBackend.from_env()
        text = backend.generate(system_prompt, user_prompt)
    """

    def __init__(self, backend_impl):
        self._impl = backend_impl

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------
    @classmethod
    def from_env(cls, env_path: str = "secrets.env") -> "LLMBackend":
        """
        Build a backend from environment variables (loaded from env_path).
        Priority: environment already set > file.
        """
        from dotenv import load_dotenv
        load_dotenv(env_path)

        backend_type = os.getenv("LLM_BACKEND", "hf_api").lower()
        model_key    = os.getenv("LLM_MODEL_KEY", "jsl_medllama_3_8b")
        hf_token     = os.getenv("HF_TOKEN", "")
        local_dir    = os.getenv("LLM_LOCAL_DIR", "/user/suvinava.basak/u19797/CoDHy_v2/models")
        device       = os.getenv("LLM_DEVICE", None)  # None → auto-detect

        if model_key not in MODEL_REGISTRY:
            raise ValueError(
                f"Unknown LLM_MODEL_KEY='{model_key}'. "
                f"Valid keys: {list(MODEL_REGISTRY.keys())}"
            )

        info = MODEL_REGISTRY[model_key]
        repo_id      = info["repo_id"]
        gated        = info["gated"]
        # Deliberately NOT reading an HF_PROVIDER env var here anymore.
        # A previous version of this code let HF_PROVIDER override the
        # registry's provider, and a stale HF_PROVIDER=auto left in
        # secrets.env from earlier setup silently broke the model:provider
        # string (e.g. sent "...:auto" instead of "...:featherless-ai"),
        # causing a real "provider or policy ... 'auto' is not valid" error.
        # The registry's provider is already verified per-model — use it
        # directly. If you genuinely need to override it, use
        # LLMBackend.from_config(provider=...) instead of an env var.
        provider     = info["provider"]
        call_style   = info["call_style"]
        max_tok_cap  = info["max_tokens_tested"]
        print(f"[LLMBackend] Selected model: {model_key} ({repo_id})")
        print(f"             Notes: {info['notes']}")
        print(f"             Backend: {backend_type}")
        if gated is None:
            print(f"             ⚠️  Gating status for this model is NOT verified — check "
                  f"https://huggingface.co/{repo_id} directly.")

        if backend_type == "ollama":
            ollama_model = os.getenv("OLLAMA_MODEL", "llama3.1")
            ollama_url   = os.getenv("OLLAMA_BASE_URL", "http://localhost:11434/api/generate")
            impl = _OllamaBackend(ollama_model, ollama_url)
            print(f"             Ollama model: {ollama_model} @ {ollama_url}")

        elif backend_type == "local":
            model_dir = os.path.join(local_dir, model_key)
            if not os.path.isdir(model_dir):
                raise FileNotFoundError(
                    f"Model directory not found: {model_dir}\n"
                    f"Run login_node_preprocess.py first to download the model."
                )
            impl = _LocalHFBackend(model_dir, device=device)
            print(f"             Local path: {model_dir}")

        elif backend_type == "hf_api":
            if gated and not hf_token:
                raise ValueError(
                    f"Model '{model_key}' ({repo_id}) is gated but no HF_TOKEN was set."
                )
            if not hf_token:
                print("             ⚠️  No HF_TOKEN set — required for hf_api, and gating "
                      "status for some models isn't verified, so this may fail even for "
                      "models not marked gated.")
            impl = _HFInferenceBackend(
                repo_id, token=hf_token, provider=provider,
                call_style=call_style, max_tokens_cap=max_tok_cap,
            )
            print(f"             HF Inference: provider={provider}, call_style={call_style}, "
                  f"max_tokens_cap={max_tok_cap} (repo: {repo_id})")

        else:
            raise ValueError(
                f"Unknown LLM_BACKEND='{backend_type}'. Use 'local', 'ollama', or 'hf_api'."
            )

        return cls(impl)

    @classmethod
    def from_config(
        cls,
        model_key: str = "jsl_medllama_3_8b",
        backend_type: str = "hf_api",
        local_dir: str = "/user/suvinava.basak/u19797/CoDHy_v2/models",
        hf_token: str = "",
        device: Optional[str] = None,
        ollama_model: str = "llama3.1",
        ollama_url: str = "http://localhost:11434/api/generate",
        provider: Optional[str] = None,
    ) -> "LLMBackend":
        """
        Programmatic alternative to from_env().
        provider : str, optional
            HF Inference Provider to use (only relevant for backend_type="hf_api").
            Defaults to the registry's provider if not given.
        """
        if model_key not in MODEL_REGISTRY:
            raise ValueError(f"Unknown model key '{model_key}'.")
        info = MODEL_REGISTRY[model_key]
        repo_id     = info["repo_id"]
        gated       = info["gated"]
        call_style  = info["call_style"]
        max_tok_cap = info["max_tokens_tested"]
        if provider is None:
            provider = info["provider"]
        print(f"[LLMBackend] {model_key} ({repo_id}) via {backend_type}")
        if gated is None:
            print(f"             ⚠️  Gating status for this model is NOT verified — check "
                  f"https://huggingface.co/{repo_id} directly.")

        if backend_type == "local":
            model_dir = os.path.join(local_dir, model_key)
            impl = _LocalHFBackend(model_dir, device=device)
        elif backend_type == "ollama":
            impl = _OllamaBackend(ollama_model, ollama_url)
        elif backend_type == "hf_api":
            impl = _HFInferenceBackend(
                repo_id, token=hf_token, provider=provider,
                call_style=call_style, max_tokens_cap=max_tok_cap,
            )
        else:
            raise ValueError(f"Unknown backend '{backend_type}'.")
        return cls(impl)

    # ------------------------------------------------------------------
    # Core API
    # ------------------------------------------------------------------
    def generate(
        self,
        system_prompt: str,
        user_prompt: str,
        temperature: float = 0.0,
        max_tokens: int = 3072,
        json_mode: bool = False,
    ) -> str:
        """Generate text; returns raw string output from the model."""
        return self._impl.generate(
            system_prompt, user_prompt,
            temperature=temperature,
            max_tokens=max_tokens,
            json_mode=json_mode,
        )

    def check_connection(self) -> bool:
        return self._impl.check_connection()

    # ------------------------------------------------------------------
    # Convenience JSON helpers (used by agents)
    # ------------------------------------------------------------------
    @staticmethod
    def extract_json_list(text: str) -> list:
        """
        Robustly extract a JSON list from messy LLM output.

        Strips a leading <think>...</think> reasoning block (DeepSeek-R1-style
        "thinking" models) BEFORE attempting bracket search. This matters
        because the old bracket-search fallback (`\\[.*\\]` with DOTALL) is
        greedy across the WHOLE text — if the reasoning block happens to
        contain any unrelated "[...]" (e.g. a bracketed pathway name), the
        greedy match spans from that stray bracket all the way to the real
        closing bracket, producing invalid JSON. Verified this failure mode
        empirically with a reasoning-block test case before writing this fix.

        Robust to a missing opening <think> tag: matches "everything up to
        and including the first </think>" rather than requiring a paired
        <think>...</think>, since a real HF discussion thread on
        deepseek-ai/DeepSeek-R1-Distill-Llama-8B reported the opening tag
        being dropped in some model revisions while the closing tag remained.
        """
        if '</think>' in text.lower():
            idx = text.lower().find('</think>') + len('</think>')
            text = text[idx:]
        text = text.strip()

        text = re.sub(r'```json', '', text)
        text = re.sub(r'```', '', text).strip()
        try:
            data = json.loads(text)
            if isinstance(data, list):
                return [x for x in data if isinstance(x, dict)]
            if isinstance(data, dict):
                for v in data.values():
                    if isinstance(v, list):
                        return [x for x in v if isinstance(x, dict)]
                return [data]
        except Exception:
            pass
        try:
            m = re.search(r'\[.*\]', text, re.DOTALL)
            if m:
                return [x for x in json.loads(m.group()) if isinstance(x, dict)]
        except Exception:
            pass
        return []

    @staticmethod
    def extract_json_dict(text: str) -> dict:
        """Robustly extract a JSON dict from messy LLM output."""
        if '</think>' in text.lower():
            idx = text.lower().find('</think>') + len('</think>')
            text = text[idx:]
        text = text.strip()

        text = re.sub(r'```json', '', text)
        text = re.sub(r'```', '', text).strip()
        try:
            data = json.loads(text)
            if isinstance(data, dict):
                return data
        except Exception:
            pass
        try:
            m = re.search(r'\{.*\}', text, re.DOTALL)
            if m:
                return json.loads(m.group())
        except Exception:
            pass
        return {}