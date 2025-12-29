# Ollama API Public Endpoint Guide

## Overview

This guide provides comprehensive instructions for accessing and using our public Ollama API endpoint for large language model inference. Our API provides access to multiple state-of-the-art language models for various use cases.

**Base URL:** `https://bmblx.bmi.osumc.edu/ollama/`

---

## 🔐 Authentication

**Currently:** No authentication required (public access)

To enable authentication in the future, the API supports:
- HTTP Basic Authentication
- Bearer Token Authentication
- API Key via headers

---

## 🤖 Available Models

| Model | Size | Parameters | Specialization | Best For |
|-------|------|------------|----------------|----------|
| **qwen3:0.6b** | 523 MB | 751M | Fast, lightweight | Quick responses, simple tasks, testing |
| **deepseek-r1:8b** | 5.2 GB | 8.2B | Reasoning-focused | Logic, problem-solving, analysis |
| **gpt-oss:20b** | 13.8 GB | 20.9B | General purpose | Balanced performance, versatile tasks |
| **qwen3:30b** | 18.6 GB | 30.5B | High-quality responses | Complex reasoning, detailed answers |
| **qwen3-coder:latest** | 18.6 GB | 30.5B | Code generation | Programming, debugging, code explanation |
| **qwen3:235b** | 142 GB | 235B | State-of-the-art | Research, complex analysis, long-form content |

### Model Selection Guide

- **Development/Testing**: Use `qwen3:0.6b` for fast iteration
- **Code Tasks**: Use `qwen3-coder:latest` for programming assistance
- **Reasoning**: Use `deepseek-r1:8b` for logical analysis
- **Production (Balanced)**: Use `qwen3:30b` for quality + speed
- **Maximum Quality**: Use `qwen3:235b` for best possible output

---

## 📡 API Endpoints

### 1. List Available Models

Get a list of all available models and their details.

**Endpoint:** `GET /api/tags`

**Example Request:**
```bash
curl https://bmblx.bmi.osumc.edu/ollama/api/tags
```

**Example Response:**
```json
{
  "models": [
    {
      "name": "qwen3:0.6b",
      "size": 522653767,
      "modified_at": "2025-10-31T14:05:58.11931-04:00",
      "details": {
        "parameter_size": "751.63M",
        "quantization_level": "Q4_K_M",
        "family": "qwen3"
      }
    }
  ]
}
```

---

### 2. Generate Text (Completion)

Generate text completion from a prompt.

**Endpoint:** `POST /api/generate`

**Headers:**
```
Content-Type: application/json
```

**Request Body:**
```json
{
  "model": "qwen3:30b",
  "prompt": "Explain quantum computing in simple terms",
  "stream": false
}
```

**Example Request (cURL):**
```bash
curl -X POST https://bmblx.bmi.osumc.edu/ollama/api/generate \
  -H "Content-Type: application/json" \
  -d '{
    "model": "qwen3:30b",
    "prompt": "What is artificial intelligence?",
    "stream": false
  }'
```

**Example Request (Python):**
```python
import requests

url = "https://bmblx.bmi.osumc.edu/ollama/api/generate"
payload = {
    "model": "qwen3:30b",
    "prompt": "Explain machine learning in simple terms",
    "stream": False
}

response = requests.post(url, json=payload)
result = response.json()
print(result['response'])
```

**Example Request (JavaScript/Node.js):**
```javascript
const response = await fetch('https://bmblx.bmi.osumc.edu/ollama/api/generate', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    model: 'qwen3:30b',
    prompt: 'What is quantum computing?',
    stream: false
  })
});

const data = await response.json();
console.log(data.response);
```

**Response:**
```json
{
  "model": "qwen3:30b",
  "created_at": "2025-10-31T18:30:00Z",
  "response": "Artificial intelligence (AI) refers to...",
  "done": true,
  "total_duration": 5000000000,
  "load_duration": 1000000000,
  "prompt_eval_count": 15,
  "eval_count": 50
}
```

---

### 3. Chat Completion

Generate conversational responses with message history.

**Endpoint:** `POST /api/chat`

**Request Body:**
```json
{
  "model": "qwen3:30b",
  "messages": [
    {
      "role": "system",
      "content": "You are a helpful AI assistant specialized in science."
    },
    {
      "role": "user",
      "content": "Explain photosynthesis"
    }
  ],
  "stream": false
}
```

**Example Request (cURL):**
```bash
curl -X POST https://bmblx.bmi.osumc.edu/ollama/api/chat \
  -H "Content-Type: application/json" \
  -d '{
    "model": "qwen3:30b",
    "messages": [
      {"role": "system", "content": "You are a helpful assistant."},
      {"role": "user", "content": "Write a haiku about technology"}
    ],
    "stream": false
  }'
```

**Example Request (Python):**
```python
import requests

url = "https://bmblx.bmi.osumc.edu/ollama/api/chat"
payload = {
    "model": "qwen3:30b",
    "messages": [
        {"role": "system", "content": "You are a helpful AI assistant."},
        {"role": "user", "content": "What is machine learning?"}
    ],
    "stream": False
}

response = requests.post(url, json=payload)
result = response.json()
print(result['message']['content'])
```

**Multi-turn Conversation Example:**
```python
import requests

url = "https://bmblx.bmi.osumc.edu/ollama/api/chat"
conversation = [
    {"role": "system", "content": "You are a helpful assistant."}
]

# First message
conversation.append({"role": "user", "content": "What is Python?"})
response = requests.post(url, json={"model": "qwen3:30b", "messages": conversation, "stream": False})
assistant_reply = response.json()['message']['content']
conversation.append({"role": "assistant", "content": assistant_reply})
print(f"Assistant: {assistant_reply}")

# Follow-up message
conversation.append({"role": "user", "content": "Can you show me an example?"})
response = requests.post(url, json={"model": "qwen3:30b", "messages": conversation, "stream": False})
assistant_reply = response.json()['message']['content']
print(f"Assistant: {assistant_reply}")
```

**Response:**
```json
{
  "model": "qwen3:30b",
  "created_at": "2025-10-31T18:30:00Z",
  "message": {
    "role": "assistant",
    "content": "Photosynthesis is the process by which plants..."
  },
  "done": true
}
```

---

### 4. Streaming Responses

Get real-time streaming responses as the model generates text.

**Endpoint:** `POST /api/generate` or `POST /api/chat`

**Set `"stream": true` in the request body**

**Example Request (cURL):**
```bash
curl -X POST https://bmblx.bmi.osumc.edu/ollama/api/generate \
  -H "Content-Type: application/json" \
  -d '{
    "model": "qwen3:30b",
    "prompt": "Write a short story about space",
    "stream": true
  }'
```

**Example Request (Python with streaming):**
```python
import requests
import json

url = "https://bmblx.bmi.osumc.edu/ollama/api/generate"
payload = {
    "model": "qwen3:30b",
    "prompt": "Write a poem about AI",
    "stream": True
}

response = requests.post(url, json=payload, stream=True)

for line in response.iter_lines():
    if line:
        chunk = json.loads(line)
        if 'response' in chunk:
            print(chunk['response'], end='', flush=True)
```

**Example Request (JavaScript):**
```javascript
const response = await fetch('https://bmblx.bmi.osumc.edu/ollama/api/generate', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    model: 'qwen3:30b',
    prompt: 'Tell me about the universe',
    stream: true
  })
});

const reader = response.body.getReader();
const decoder = new TextDecoder();

while (true) {
  const { done, value } = await reader.read();
  if (done) break;
  
  const chunk = decoder.decode(value);
  const lines = chunk.split('\n').filter(line => line.trim());
  
  for (const line of lines) {
    const data = JSON.parse(line);
    if (data.response) {
      process.stdout.write(data.response);
    }
  }
}
```

**Streaming Response Format:**
Each line is a separate JSON object:
```json
{"model":"qwen3:30b","response":"Once","done":false}
{"model":"qwen3:30b","response":" upon","done":false}
{"model":"qwen3:30b","response":" a","done":false}
...
{"model":"qwen3:30b","response":"","done":true}
```

---

### 5. Check Running Models

See which models are currently loaded in memory.

**Endpoint:** `GET /api/ps`

**Example Request:**
```bash
curl https://bmblx.bmi.osumc.edu/ollama/api/ps
```

**Example Response:**
```json
{
  "models": [
    {
      "name": "qwen3:30b",
      "size": 18556699314,
      "expires_at": "2318-02-10T14:38:08.312955441-05:00",
      "size_vram": 18556699314
    }
  ]
}
```

---

## 🎛️ Advanced Parameters

### Common Parameters for `/api/generate`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model` | string | required | Model name to use |
| `prompt` | string | required | The prompt to generate from |
| `stream` | boolean | false | Enable streaming responses |
| `temperature` | float | 0.8 | Randomness (0.0-2.0). Lower = more focused |
| `top_p` | float | 0.9 | Nucleus sampling threshold |
| `top_k` | integer | 40 | Top-k sampling |
| `num_predict` | integer | -1 | Max tokens to generate (-1 = unlimited) |
| `stop` | array | null | Stop sequences |
| `seed` | integer | random | Random seed for reproducibility |
| `repeat_penalty` | float | 1.1 | Penalize repetition |
| `presence_penalty` | float | 0.0 | Penalize tokens based on presence |
| `frequency_penalty` | float | 0.0 | Penalize tokens based on frequency |

**Example with Advanced Parameters:**
```json
{
  "model": "qwen3:30b",
  "prompt": "Explain neural networks",
  "stream": false,
  "temperature": 0.7,
  "top_p": 0.9,
  "top_k": 50,
  "num_predict": 500,
  "repeat_penalty": 1.2,
  "stop": ["\n\n", "###"]
}
```

### Temperature Guide

- **0.0-0.3**: Very focused, deterministic (good for factual questions, code)
- **0.4-0.7**: Balanced creativity and coherence (general use)
- **0.8-1.2**: More creative and varied (creative writing, brainstorming)
- **1.3-2.0**: Highly creative, potentially chaotic (experimental)

---

## 💻 Code Examples by Use Case

### Use Case 1: Code Generation

```python
import requests

def generate_code(description, language="Python"):
    url = "https://bmblx.bmi.osumc.edu/ollama/api/generate"
    payload = {
        "model": "qwen3-coder:latest",
        "prompt": f"Write {language} code for: {description}\n\nProvide only the code without explanation.",
        "stream": False,
        "temperature": 0.3
    }
    response = requests.post(url, json=payload)
    return response.json()['response']

# Example
code = generate_code("a function that calculates fibonacci numbers")
print(code)
```

### Use Case 2: Text Summarization

```python
def summarize_text(text, max_sentences=3):
    url = "https://bmblx.bmi.osumc.edu/ollama/api/generate"
    payload = {
        "model": "qwen3:30b",
        "prompt": f"Summarize the following text in {max_sentences} sentences:\n\n{text}",
        "stream": False,
        "temperature": 0.5
    }
    response = requests.post(url, json=payload)
    return response.json()['response']
```

### Use Case 3: Question Answering with Context

```python
def qa_with_context(context, question):
    url = "https://bmblx.bmi.osumc.edu/ollama/api/chat"
    payload = {
        "model": "qwen3:30b",
        "messages": [
            {"role": "system", "content": f"Answer questions based on this context: {context}"},
            {"role": "user", "content": question}
        ],
        "stream": False
    }
    response = requests.post(url, json=payload)
    return response.json()['message']['content']
```

### Use Case 4: Translation

```python
def translate(text, target_language):
    url = "https://bmblx.bmi.osumc.edu/ollama/api/generate"
    payload = {
        "model": "qwen3:30b",
        "prompt": f"Translate the following text to {target_language}:\n\n{text}",
        "stream": False
    }
    response = requests.post(url, json=payload)
    return response.json()['response']
```

### Use Case 5: Reasoning and Analysis

```python
def analyze_problem(problem):
    url = "https://bmblx.bmi.osumc.edu/ollama/api/chat"
    payload = {
        "model": "deepseek-r1:8b",
        "messages": [
            {"role": "system", "content": "You are a logical reasoning assistant. Break down problems step by step."},
            {"role": "user", "content": problem}
        ],
        "stream": False
    }
    response = requests.post(url, json=payload)
    return response.json()['message']['content']
```

---

## 🔧 Integration with Popular Tools

### LangChain (Python)

```python
from langchain.llms import Ollama

llm = Ollama(
    base_url="https://bmblx.bmi.osumc.edu/ollama",
    model="qwen3:30b"
)

response = llm("What is the theory of relativity?")
print(response)
```

### OpenAI Python Client (Compatible)

```python
from openai import OpenAI

client = OpenAI(
    base_url="https://bmblx.bmi.osumc.edu/ollama/v1",
    api_key="ollama"  # Required but ignored
)

response = client.chat.completions.create(
    model="qwen3:30b",
    messages=[
        {"role": "user", "content": "Explain photosynthesis"}
    ]
)

print(response.choices[0].message.content)
```

### LiteLLM (Multi-provider wrapper)

```python
from litellm import completion

response = completion(
    model="ollama/qwen3:30b",
    messages=[{"role": "user", "content": "Hello, how are you?"}],
    api_base="https://bmblx.bmi.osumc.edu/ollama"
)
print(response.choices[0].message.content)
```

### Continue.dev (VS Code Extension)

Add to your `config.json`:
```json
{
  "models": [
    {
      "title": "Qwen3 30B",
      "provider": "ollama",
      "model": "qwen3:30b",
      "apiBase": "https://bmblx.bmi.osumc.edu/ollama"
    },
    {
      "title": "Qwen3 Coder",
      "provider": "ollama",
      "model": "qwen3-coder:latest",
      "apiBase": "https://bmblx.bmi.osumc.edu/ollama"
    }
  ]
}
```

### Cursor IDE

In Cursor settings:
```json
{
  "cursor.aiProvider": "ollama",
  "cursor.ollamaBaseUrl": "https://bmblx.bmi.osumc.edu/ollama",
  "cursor.ollamaModel": "qwen3-coder:latest"
}
```

---

## 🧪 Testing with API Clients

### Postman

1. Create a new POST request
2. URL: `https://bmblx.bmi.osumc.edu/ollama/api/generate`
3. Headers: `Content-Type: application/json`
4. Body (raw JSON):
```json
{
  "model": "qwen3:30b",
  "prompt": "Your prompt here",
  "stream": false
}
```
5. Click Send

### Insomnia

1. Create a new POST request
2. URL: `https://bmblx.bmi.osumc.edu/ollama/api/generate`
3. Headers: `Content-Type: application/json`
4. Body: Select JSON
```json
{
  "model": "qwen3:30b",
  "prompt": "Your prompt here",
  "stream": false
}
```
5. Click Send

### Thunder Client (VS Code)

1. New Request → POST
2. URL: `https://bmblx.bmi.osumc.edu/ollama/api/generate`
3. Headers: `Content-Type: application/json`
4. Body: JSON
```json
{
  "model": "qwen3:30b",
  "prompt": "Your prompt here",
  "stream": false
}
```

---

## ⚡ Performance Characteristics

### Expected Response Times

| Model | Tokens/Second | First Token Latency | Best For |
|-------|---------------|---------------------|----------|
| qwen3:0.6b | ~50-100 | <1s | Real-time applications |
| deepseek-r1:8b | ~30-50 | 1-2s | Balanced performance |
| qwen3:30b | ~20-40 | 2-3s | Production workloads |
| qwen3-coder:latest | ~20-40 | 2-3s | Code generation |
| qwen3:235b | ~30-40 | 3-5s | Maximum quality |

**Note:** Response times vary based on:
- Prompt length
- Generated token count
- Server load
- Model warmup state

### Token Estimation

- **Average**: ~4 characters = 1 token
- **English words**: ~1.3 tokens per word
- **Code**: ~1.5 tokens per word
- **Estimation formula**: `tokens ≈ characters / 4`

---

## 🛡️ Best Practices

### 1. Model Selection

- Use **qwen3:0.6b** for prototyping and testing
- Use **qwen3-coder:latest** for all code-related tasks
- Use **deepseek-r1:8b** for logical reasoning and problem-solving
- Use **qwen3:30b** for production (best balance of quality and speed)
- Use **qwen3:235b** only when absolute best quality is needed

### 2. Prompt Engineering

**Be Clear and Specific:**
```json
// Bad
{"prompt": "write code"}

// Good
{"prompt": "Write a Python function that takes a list of integers and returns the sum of all even numbers. Include error handling for non-integer inputs."}
```

**Use System Messages in Chat:**
```json
{
  "messages": [
    {"role": "system", "content": "You are an expert Python developer. Provide clean, well-documented code with type hints."},
    {"role": "user", "content": "Create a function to validate email addresses"}
  ]
}
```

**Set Appropriate Temperature:**
- Creative tasks (writing, brainstorming): 0.7-1.0
- Balanced tasks (general Q&A): 0.5-0.7
- Factual tasks (code, math): 0.1-0.3

### 3. Optimize Token Usage

```python
# Limit output length
{
  "num_predict": 500,  # Limit to 500 tokens
  "stop": ["\n\n", "###"]  # Stop at double newline or ###
}
```

### 4. Error Handling

```python
import requests
import time

def api_call_with_retry(url, payload, max_retries=3):
    for attempt in range(max_retries):
        try:
            response = requests.post(url, json=payload, timeout=60)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.Timeout:
            if attempt == max_retries - 1:
                raise
            time.sleep(2 ** attempt)  # Exponential backoff
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")
            if attempt == max_retries - 1:
                raise
            time.sleep(2 ** attempt)
```

### 5. Use Streaming for Better UX

For long responses, always use streaming to show progress:

```python
import requests
import json

def stream_response(prompt, model="qwen3:30b"):
    url = "https://bmblx.bmi.osumc.edu/ollama/api/generate"
    payload = {"model": model, "prompt": prompt, "stream": True}
    
    response = requests.post(url, json=payload, stream=True)
    
    full_response = ""
    for line in response.iter_lines():
        if line:
            chunk = json.loads(line)
            if 'response' in chunk:
                print(chunk['response'], end='', flush=True)
                full_response += chunk['response']
    
    return full_response
```

---

## 🚨 Error Handling

### Common HTTP Status Codes

| Code | Meaning | Solution |
|------|---------|----------|
| 200 | Success | Request completed successfully |
| 400 | Bad Request | Check JSON format and parameters |
| 404 | Not Found | Model doesn't exist, verify model name |
| 500 | Server Error | Server issue, try again later or use smaller model |
| 503 | Service Unavailable | Model loading or server overloaded, retry with backoff |

### Example Error Responses

**Model Not Found:**
```json
{
  "error": "model 'qwen3:999b' not found, try pulling it first"
}
```

**Invalid Request:**
```json
{
  "error": "invalid request: missing required field 'model'"
}
```

---

## 📊 Rate Limits and Usage Guidelines

### Current Limits
- **No hard rate limits** currently enforced
- Please be considerate of shared resources
- Recommended: Max 10 concurrent requests per user

### Fair Usage Guidelines

1. **Use appropriate models** - Don't use qwen3:235b for simple tasks
2. **Set reasonable token limits** - Use `num_predict` to limit output
3. **Implement backoff** - If getting 503 errors, wait before retrying
4. **Cache responses** - Don't repeatedly request identical prompts
5. **Use streaming** - For better resource utilization

---

## 🔍 Troubleshooting

### Issue: Slow Response Times

**Solutions:**
- Use a smaller model (qwen3:0.6b or qwen3:30b)
- Reduce `num_predict` value
- Check if server is under heavy load
- Use streaming for better perceived performance

### Issue: Empty or Incomplete Responses

**Solutions:**
- Check if `num_predict` is too low
- Verify prompt is clear and specific
- Try increasing `temperature` slightly
- Ensure model is not being unloaded (check `/api/ps`)

### Issue: Repetitive Output

**Solutions:**
- Increase `repeat_penalty` to 1.2-1.5
- Lower `temperature` to 0.3-0.5
- Add stop sequences
- Rephrase prompt to be more specific

### Issue: Model Not Available

**Solutions:**
- Check `/api/tags` for available models
- Verify model name spelling
- Wait a few seconds for model to load
- Try a different model

---

## 📞 Support and Resources

### Getting Help

- **API Status**: Check `/api/ps` for loaded models
- **Model List**: Check `/api/tags` for available models
- **Documentation**: This README

### Reporting Issues

When reporting issues, please include:
1. Model name and version
2. Full request (JSON payload)
3. Response or error message
4. Timestamp of the issue

---

## 📈 Changelog

### Version 1.1 (2025-10-31)
- Added 5 new models: qwen3:30b, qwen3-coder, deepseek-r1:8b, gpt-oss:20b
- Updated performance characteristics
- Added advanced parameter documentation
- Expanded code examples

### Version 1.0 (2025-10-31)
- Initial public release
- Models: qwen3:0.6b, qwen3:235b
- Endpoints: /api/tags, /api/generate, /api/chat, /api/ps

---

## 📜 Terms of Use

- This API is provided for **research and educational purposes**
- Please use responsibly and ethically
- Do not abuse the service with excessive requests
- Content generated by AI should be reviewed before use
- No warranty or guarantee of availability

---

## 🌟 Quick Start Examples

### Python Quick Start

```python
import requests

# Simple text generation
response = requests.post(
    "https://bmblx.bmi.osumc.edu/ollama/api/generate",
    json={
        "model": "qwen3:30b",
        "prompt": "What is machine learning?",
        "stream": False
    }
)
print(response.json()['response'])
```

### JavaScript Quick Start

```javascript
// Simple text generation
const response = await fetch('https://bmblx.bmi.osumc.edu/ollama/api/generate', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    model: 'qwen3:30b',
    prompt: 'What is machine learning?',
    stream: false
  })
});

const data = await response.json();
console.log(data.response);
```

### cURL Quick Start

```bash
# Simple text generation
curl -X POST https://bmblx.bmi.osumc.edu/ollama/api/generate \
  -H "Content-Type: application/json" \
  -d '{
    "model": "qwen3:30b",
    "prompt": "What is machine learning?",
    "stream": false
  }'
```

---

**Last Updated:** October 31, 2025  
**API Version:** 1.1  
**Base URL:** `https://bmblx.bmi.osumc.edu/ollama/`  
**Available Models:** 6 (0.6B to 235B parameters)
