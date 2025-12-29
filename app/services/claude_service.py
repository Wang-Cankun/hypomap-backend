"""
Claude AI Service for HypoMap - biological analysis and hypothesis generation
"""
import os
import anthropic
from typing import Optional


class ClaudeService:
    """Service for interacting with Claude AI for biological analysis"""

    def __init__(self):
        self.api_key = os.getenv("ANTHROPIC_API_KEY", "")
        self.client = None
        if self.api_key:
            self.client = anthropic.Anthropic(api_key=self.api_key)

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

    async def analyze_cluster(self, context: str, question: str) -> str:
        """
        Analyze a cluster using Claude AI

        Args:
            context: Formatted string with DEG, regulon, and CCC data
            question: User's question about the cluster

        Returns:
            AI-generated analysis response
        """
        if not self.client:
            return "AI service not configured. Please set ANTHROPIC_API_KEY environment variable."

        try:
            message = self.client.messages.create(
                model="claude-sonnet-4-20250514",
                max_tokens=2048,
                system=self.system_prompt,
                messages=[
                    {
                        "role": "user",
                        "content": f"""Here is the analysis context for a single-cell RNA-seq cluster:

{context}

---

Question: {question}

Please provide a detailed, scientifically accurate response based on the data provided above."""
                    }
                ]
            )

            return message.content[0].text

        except anthropic.APIConnectionError:
            return "Error: Unable to connect to Claude API. Please check your internet connection."
        except anthropic.RateLimitError:
            return "Error: API rate limit exceeded. Please try again later."
        except anthropic.APIStatusError as e:
            return f"Error: API returned status {e.status_code}: {e.message}"
        except Exception as e:
            return f"Error generating response: {str(e)}"

    async def generate_hypothesis(self, context: str) -> str:
        """Generate biological hypotheses from cluster data"""
        question = """Based on the gene expression patterns, transcription factors, and cell-cell communication data provided:

1. What is the most likely biological identity/function of this cell cluster?
2. Generate 3-5 testable hypotheses about this cluster's role in the tissue
3. What validation experiments would you recommend?
4. Are there any disease implications based on the marker genes?

Please be specific and cite the genes/pathways that support each hypothesis."""

        return await self.analyze_cluster(context, question)

    async def identify_cell_type(self, context: str) -> str:
        """Identify cell type from marker genes"""
        question = """Based on the differentially expressed genes provided:

1. What cell type is this cluster most likely to represent?
2. List the key marker genes that support this identification
3. Provide a confidence level (high/medium/low) and explain why
4. Are there any alternative cell type interpretations?
5. What additional markers would help confirm the identity?"""

        return await self.analyze_cluster(context, question)

    async def analyze_pathways(self, context: str) -> str:
        """Analyze active pathways in a cluster"""
        question = """Analyze the biological pathways active in this cluster:

1. What are the main signaling pathways represented by the DEGs?
2. What metabolic pathways appear to be active?
3. Are there any developmental or differentiation pathways?
4. How do the transcription factors relate to these pathways?
5. What functional state does this suggest for the cells?"""

        return await self.analyze_cluster(context, question)
