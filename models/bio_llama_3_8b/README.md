---
license: other
library_name: transformers
tags:
- reasoning
- chain-of-thought
- medical
- Healthcare & Lifesciences
- BioMed
base_model: ContactDoctor/Bio-Medical-Llama-3-8B
thumbnail: https://collaiborate.com/logo/logo-blue-bg-1.png
model-index:
- name: Bio-Medical-Llama-3-8B-CoT-012025
  results: []
datasets:
- collaiborateorg/BioMedData-CoT
---

# Bio-Medical-Llama-3-8B-CoT-012025

![image/jpeg](https://cdn-uploads.huggingface.co/production/uploads/653f5b93cd52f288490edc83/zPMUugzfOiwTiRw88jm7T.jpeg)

This model, **Bio-Medical-Llama-3-8B-CoT-012025**, is a fine-tuned extension of the original [Bio-Medical-Llama-3-8B](https://huggingface.co/ContactDoctor/Bio-Medical-Llama-3-8B) and [Deepseek's Distilled Llama 8B](https://huggingface.co/deepseek-ai/DeepSeek-R1-Distill-Llama-8B) version, now equipped with advanced reasoning capabilities using a Chain-of-Thought (COT) instruction set. This enhancement builds upon our commitment to delivering state-of-the-art, specialized LLMs for the healthcare and life sciences domains.

## Model Details

Model Name: Bio-Medical-Llama-3-8B-CoT-012025

Base Model: Bio-Medical-Llama-3-8B

Parameter Count: 8 billion

Training Data: Extended dataset comprising high-quality biomedical data with a focus on reasoning-intensive tasks.

Number of Entries in Original Dataset: 600K+, Extension Dataset: 25K+

Dataset Composition: The dataset integrates diverse and reasoning-centric biomedical queries and tasks, ensuring robust Chain-of-Thought performance. It includes both synthetic and manually curated examples tailored to clinical, diagnostic, and research-oriented scenarios.

## Model Description

**Bio-Medical-Llama-3-8B-CoT-012025** represents a leap forward in AI-driven reasoning for the healthcare and life sciences sectors. By incorporating Chain-of-Thought fine-tuning, the model excels at handling complex, multi-step reasoning tasks, making it ideal for scenarios requiring critical thinking and nuanced understanding.

Key Features:
- **Enhanced Reasoning Abilities**: Trained specifically to perform multi-step reasoning and provide accurate, contextually rich responses.
- **Compact Model Sizes for Versatility**: Includes 1B, 3B, and 8B variants optimized for edge devices and high-performance systems alike.
- **Specialized Training Focus**: Developed using datasets designed to address the unique challenges of biomedical reasoning and problem-solving.

## Evaluation Metrics

**Bio-Medical-Llama-3-8B-CoT-012025** demonstrates significant advancements in multi-step reasoning tasks, surpassing its predecessor on HLS benchmarks. However, it is important to note that its performance on evaluation tasks involving multiple-choice questions may be less robust, as it is specifically optimized for reasoning-based challenges.

## Intended Uses & Limitations

Bio-Medical-Llama-3-8B-CoT-012025 is designed for applications requiring high levels of reasoning within the biomedical field, including:

1. **Clinical Reasoning**: Supporting healthcare professionals in diagnostic and treatment planning.
2. **Medical Research**: Assisting in hypothesis generation, literature synthesis, and data interpretation.
3. **Educational Tools**: Providing medical students and professionals with advanced training simulations and problem-solving support.

### Limitations and Ethical Considerations

> **Biases**: While efforts were made to minimize bias during training, some biases inherent in the training data may persist.

> **Accuracy**: This model’s reasoning is based on training data and may not always be up-to-date or contextually perfect. Users should verify critical outputs against authoritative sources.

> **Ethical Use**: The model is not a substitute for professional medical judgment and should be used responsibly, particularly in clinical decision-making.

## How to Use

```python
import transformers
import torch

model_id = "ContactDoctor/Bio-Medical-Llama-3-8B-CoT-012025"

pipeline = transformers.pipeline(
    "text-generation",
    model=model_id,
    model_kwargs={"torch_dtype": torch.bfloat16},
    device_map="auto",
)

messages = [
    {"role": "system", "content": "You are an expert trained on healthcare and biomedical reasoning."},
    {"role": "user", "content": "What are the differential diagnoses for a 45-year-old male presenting with chest pain?"},
]

prompt = pipeline.tokenizer.apply_chat_template(
        messages,
        tokenize=False,
        add_generation_prompt=True
)

terminators = [
    pipeline.tokenizer.eos_token_id,
    pipeline.tokenizer.convert_tokens_to_ids("<|eot_id|>")
]

outputs = pipeline(
    prompt,
    max_new_tokens=256,
    eos_token_id=terminators,
    do_sample=True,
    temperature=0.6,
    top_p=0.9,
)
print(outputs[0]["generated_text"][len(prompt):])
```

## License

This model is licensed under the [Bio-Medical-Llama-3-8B-CoT-012025 (Non-Commercial Use Only)](./LICENSE). Please review the terms and conditions before use.

### Contact Information

For further information, inquiries, or issues related to Bio-Medical-Llama-3-8B-CoT-012025, please contact:

Email: info@contactdoctor.in

Website: https://www.contactdoctor.in

### Training Hyperparameters

The following hyperparameters were used during training:
- learning_rate: 0.00015
- train_batch_size: 16
- eval_batch_size: 8
- seed: 42
- gradient_accumulation_steps: 4
- total_train_batch_size: 64
- optimizer: AdamW with betas=(0.9, 0.98) and epsilon=1e-08
- lr_scheduler_type: linear
- lr_scheduler_warmup_ratio: 0.05
- training_steps: 2500
- mixed_precision_training: Native AMP

### Framework Versions

- PEFT 0.12.0
- Transformers 4.41.0
- Pytorch 2.1.2
- Datasets 2.21.0
- Tokenizers 0.22.0

### Citation

If you use Bio-Medical-Llama-3-8B-CoT-012025 in your research or applications, please cite it as follows:

@misc{ContactDoctor_Bio-Medical-Llama-3-8B-CoT,
  author = ContactDoctor,
  title = {Bio-Medical-CoT: Advanced Reasoning for Healthcare Applications},
  year = {2025},
  howpublished = {https://huggingface.co/ContactDoctor/Bio-Medical-Llama-3-8B-CoT-012025},
}

