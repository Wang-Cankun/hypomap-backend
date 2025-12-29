"""
Ollama AI Service for HypoMap - biological analysis and hypothesis generation
Uses locally hosted Ollama at bmblx.bmi.osumc.edu
"""
import httpx
import json
from typing import Optional, AsyncGenerator

# Ollama server configuration
OLLAMA_BASE_URL = "https://bmblx.bmi.osumc.edu/ollama"

# Available models
AVAILABLE_MODELS = {
    "qwen3:30b": {
        "name": "Qwen3 30B",
        "description": "High-quality responses, best for complex biology questions",
        "size": "18.6 GB"
    },
    "gpt-oss:20b": {
        "name": "GPT-OSS 20B",
        "description": "General purpose, versatile tasks",
        "size": "13.8 GB"
    },
    "deepseek-r1:8b": {
        "name": "DeepSeek R1 8B",
        "description": "Reasoning-focused, good for analysis",
        "size": "5.2 GB"
    },
    "qwen3:0.6b": {
        "name": "Qwen3 0.6B",
        "description": "Fast, lightweight for quick responses",
        "size": "523 MB"
    }
}

DEFAULT_MODEL = "qwen3:30b"


class OllamaService:
    """Service for interacting with Ollama AI for biological analysis"""

    def __init__(self):
        self.base_url = OLLAMA_BASE_URL
        self.system_prompt = """You are an expert single-cell biologist and bioinformatician helping researchers analyze single-cell RNA sequencing data.

Your role is to:
1. Interpret gene expression patterns and marker genes
2. Identify cell types based on canonical markers
3. Explain biological pathways and their relevance
4. Analyze transcription factor regulatory networks
5. Interpret cell-cell communication patterns
6. Generate testable biological hypotheses
7. Suggest validation experiments

When analyzing data:
- Be specific about gene functions and known biology
- Reference canonical markers for cell type identification
- Consider the tissue/organ context when available
- Highlight unexpected or interesting findings
- Provide confidence levels for cell type assignments
- Suggest follow-up analyses when appropriate

Format your responses with clear sections and bullet points for readability.
Use markdown formatting for headers and emphasis."""

    async def get_available_models(self) -> dict:
        """Get list of available models"""
        return AVAILABLE_MODELS

    async def check_status(self) -> dict:
        """Check if Ollama service is available"""
        try:
            async with httpx.AsyncClient(timeout=10.0) as client:
                response = await client.get(f"{self.base_url}/api/tags")
                if response.status_code == 200:
                    data = response.json()
                    return {
                        "available": True,
                        "models": [m["name"] for m in data.get("models", [])],
                        "message": "Ollama service is ready"
                    }
        except Exception as e:
            return {
                "available": False,
                "models": [],
                "message": f"Cannot connect to Ollama: {str(e)}"
            }
        return {
            "available": False,
            "models": [],
            "message": "Ollama service unavailable"
        }

    async def analyze_cluster(
        self,
        context: str,
        question: str,
        model: str = DEFAULT_MODEL,
        stream: bool = False
    ) -> str:
        """
        Analyze a cluster using Ollama AI

        Args:
            context: Formatted string with DEG, regulon, and CCC data
            question: User's question about the cluster
            model: Model to use (default: qwen3:30b)
            stream: Whether to stream the response

        Returns:
            AI-generated analysis response
        """
        messages = [
            {
                "role": "system",
                "content": self.system_prompt
            },
            {
                "role": "user",
                "content": f"""Here is the analysis context for a single-cell RNA-seq cluster:

{context}

---

Question: {question}

Please provide a detailed, scientifically accurate response based on the data provided above."""
            }
        ]

        try:
            async with httpx.AsyncClient(timeout=120.0) as client:
                response = await client.post(
                    f"{self.base_url}/api/chat",
                    json={
                        "model": model,
                        "messages": messages,
                        "stream": False,
                        "options": {
                            "temperature": 0.7,
                            "top_p": 0.9,
                            "num_predict": 2048
                        }
                    }
                )

                if response.status_code != 200:
                    return f"Error: Ollama returned status {response.status_code}"

                data = response.json()
                return data.get("message", {}).get("content", "No response generated")

        except httpx.TimeoutException:
            return "Error: Request timed out. The model may be loading or the question is complex. Please try again."
        except httpx.ConnectError:
            return "Error: Cannot connect to Ollama server. Please check if the service is running."
        except Exception as e:
            return f"Error generating response: {str(e)}"

    async def analyze_cluster_stream(
        self,
        context: str,
        question: str,
        model: str = DEFAULT_MODEL
    ) -> AsyncGenerator[str, None]:
        """
        Stream analyze a cluster using Ollama AI

        Yields:
            Chunks of the AI-generated response
        """
        messages = [
            {
                "role": "system",
                "content": self.system_prompt
            },
            {
                "role": "user",
                "content": f"""Here is the analysis context for a single-cell RNA-seq cluster:

{context}

---

Question: {question}

Please provide a detailed, scientifically accurate response based on the data provided above."""
            }
        ]

        try:
            async with httpx.AsyncClient(timeout=120.0) as client:
                async with client.stream(
                    "POST",
                    f"{self.base_url}/api/chat",
                    json={
                        "model": model,
                        "messages": messages,
                        "stream": True,
                        "options": {
                            "temperature": 0.7,
                            "top_p": 0.9,
                            "num_predict": 2048
                        }
                    }
                ) as response:
                    async for line in response.aiter_lines():
                        if line:
                            try:
                                data = json.loads(line)
                                if "message" in data and "content" in data["message"]:
                                    yield data["message"]["content"]
                            except json.JSONDecodeError:
                                continue

        except Exception as e:
            yield f"Error: {str(e)}"

    async def chat_with_history(
        self,
        messages: list,
        model: str = DEFAULT_MODEL
    ) -> str:
        """
        Chat with conversation history

        Args:
            messages: List of message dicts with 'role' and 'content'
            model: Model to use

        Returns:
            AI response
        """
        # Prepend system message
        full_messages = [
            {"role": "system", "content": self.system_prompt}
        ] + messages

        try:
            async with httpx.AsyncClient(timeout=120.0) as client:
                response = await client.post(
                    f"{self.base_url}/api/chat",
                    json={
                        "model": model,
                        "messages": full_messages,
                        "stream": False,
                        "options": {
                            "temperature": 0.7,
                            "top_p": 0.9,
                            "num_predict": 2048
                        }
                    }
                )

                if response.status_code != 200:
                    return f"Error: Ollama returned status {response.status_code}"

                data = response.json()
                return data.get("message", {}).get("content", "No response generated")

        except Exception as e:
            return f"Error: {str(e)}"
