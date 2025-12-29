"""
AI endpoints for HypoMap - Ollama-powered biological analysis Q&A
"""
from fastapi import APIRouter, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from typing import List, Optional
from app.services.ollama_service import OllamaService, AVAILABLE_MODELS, DEFAULT_MODEL
from app.services.h5ad_service import H5ADService

router = APIRouter()

# Initialize services
ollama_service = OllamaService()
h5ad_service = H5ADService()


class AskRequest(BaseModel):
    dataset_id: str
    cluster: str
    question: str
    model: str = DEFAULT_MODEL
    include_context: List[str] = ["deg", "regulon", "ccc"]
    stream: bool = False


class ChatMessage(BaseModel):
    role: str  # "user" or "assistant"
    content: str


class ChatRequest(BaseModel):
    dataset_id: str
    cluster: str
    messages: List[ChatMessage]
    model: str = DEFAULT_MODEL
    include_context: List[str] = ["deg", "regulon", "ccc"]


class AskResponse(BaseModel):
    question: str
    answer: str
    cluster: str
    dataset_id: str
    model: str
    context_used: List[str]


@router.get("/ai/status")
async def get_ai_status():
    """Check if AI service is available"""
    status = await ollama_service.check_status()
    return status


@router.get("/ai/models")
async def get_available_models():
    """Get list of available AI models"""
    return {
        "models": AVAILABLE_MODELS,
        "default": DEFAULT_MODEL
    }


@router.post("/ai/ask")
async def ask_question(request: AskRequest):
    """Ask a question about a cluster using AI"""
    # Validate model
    if request.model not in AVAILABLE_MODELS:
        raise HTTPException(
            status_code=400,
            detail=f"Invalid model. Available: {list(AVAILABLE_MODELS.keys())}"
        )

    # Build context from available data
    context = build_cluster_context(
        request.dataset_id,
        request.cluster,
        request.include_context
    )

    # Handle streaming response
    if request.stream:
        async def generate():
            async for chunk in ollama_service.analyze_cluster_stream(
                context,
                request.question,
                model=request.model
            ):
                yield chunk

        return StreamingResponse(
            generate(),
            media_type="text/plain"
        )

    # Non-streaming response
    try:
        answer = await ollama_service.analyze_cluster(
            context,
            request.question,
            model=request.model
        )
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"AI service error: {str(e)}"
        )

    return AskResponse(
        question=request.question,
        answer=answer,
        cluster=request.cluster,
        dataset_id=request.dataset_id,
        model=request.model,
        context_used=request.include_context
    )


@router.post("/ai/chat")
async def chat_with_context(request: ChatRequest):
    """Chat with AI including conversation history"""
    # Validate model
    if request.model not in AVAILABLE_MODELS:
        raise HTTPException(
            status_code=400,
            detail=f"Invalid model. Available: {list(AVAILABLE_MODELS.keys())}"
        )

    # Build context
    context = build_cluster_context(
        request.dataset_id,
        request.cluster,
        request.include_context
    )

    # Inject context into first user message
    messages = []
    context_injected = False

    for msg in request.messages:
        if msg.role == "user" and not context_injected:
            # Add context to first user message
            messages.append({
                "role": "user",
                "content": f"""Analysis Context:
{context}

---

{msg.content}"""
            })
            context_injected = True
        else:
            messages.append({"role": msg.role, "content": msg.content})

    try:
        answer = await ollama_service.chat_with_history(
            messages,
            model=request.model
        )
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"AI service error: {str(e)}"
        )

    return {
        "response": answer,
        "model": request.model
    }


@router.get("/ai/suggested-questions/{dataset_id}/{cluster}")
async def get_suggested_questions(dataset_id: str, cluster: str):
    """Get suggested questions for a cluster"""
    questions = [
        {
            "id": "gene_function",
            "question": f"What is the function of the key marker genes in cluster {cluster}?",
            "category": "Gene Analysis"
        },
        {
            "id": "top_gene",
            "question": f"What is the highest expressed gene in cluster {cluster} and what does it indicate?",
            "category": "Gene Analysis"
        },
        {
            "id": "cell_type",
            "question": f"What cell type is cluster {cluster} likely to be based on its marker genes?",
            "category": "Cell Type Identification"
        },
        {
            "id": "pathways",
            "question": f"What biological pathways are most active in cluster {cluster} for validation?",
            "category": "Pathway Analysis"
        },
        {
            "id": "tf_regulation",
            "question": f"Which transcription factors are mainly influencing cluster {cluster}?",
            "category": "Regulatory Analysis"
        },
        {
            "id": "ccc_signals",
            "question": f"What cell-cell communication signals is cluster {cluster} sending and receiving?",
            "category": "Cell Communication"
        },
        {
            "id": "hypothesis",
            "question": f"Generate testable hypotheses about the biological role of cluster {cluster}.",
            "category": "Hypothesis Generation"
        },
        {
            "id": "validation",
            "question": f"What experiments would you recommend to validate the identity of cluster {cluster}?",
            "category": "Experimental Design"
        }
    ]
    return {"questions": questions, "cluster": cluster, "dataset_id": dataset_id}


@router.get("/ai/context/{dataset_id}/{cluster}")
async def get_cluster_context(dataset_id: str, cluster: str):
    """Get the context that would be sent to AI for a cluster"""
    context = build_cluster_context(dataset_id, cluster, ["deg", "regulon", "ccc"])
    return {"context": context, "cluster": cluster, "dataset_id": dataset_id}


def build_cluster_context(dataset_id: str, cluster: str, include: List[str]) -> str:
    """Build context string from available data sources"""
    sections = []

    # Dataset info
    try:
        info = h5ad_service.get_dataset_info(dataset_id)
        sections.append(f"""## Dataset Information
- Dataset ID: {dataset_id}
- Total cells: {info.get('n_cells', 'Unknown')}
- Total genes: {info.get('n_genes', 'Unknown')}
- Cluster being analyzed: {cluster}""")
    except Exception:
        sections.append(f"## Dataset: {dataset_id}, Cluster: {cluster}")

    # DEG data
    if "deg" in include:
        try:
            deg_data = h5ad_service.get_precomputed_deg(dataset_id, cluster, top_n=20)
            if deg_data:
                deg_text = "\n".join([
                    f"- {g['gene']}: logFC={g.get('logfoldchanges', 'N/A'):.2f}, pval={g.get('pvals_adj', 'N/A'):.2e}"
                    for g in deg_data[:20]
                ])
                sections.append(f"""## Top Differentially Expressed Genes (Cluster {cluster})
{deg_text}""")
        except Exception:
            pass

    # Regulon data
    if "regulon" in include:
        try:
            regulon_data = h5ad_service.get_regulon_data(dataset_id, cluster)
            if regulon_data:
                # Get top TFs by target count
                tf_counts = {}
                for edge in regulon_data.get('edges', []):
                    tf = edge.get('tf', edge.get('TF'))
                    tf_counts[tf] = tf_counts.get(tf, 0) + 1

                top_tfs = sorted(tf_counts.items(), key=lambda x: x[1], reverse=True)[:10]
                tf_text = "\n".join([f"- {tf}: {count} target genes" for tf, count in top_tfs])
                sections.append(f"""## Active Transcription Factors (Cluster {cluster})
{tf_text}""")
        except Exception:
            pass

    # CCC data
    if "ccc" in include:
        try:
            ccc_data = h5ad_service.get_ccc_data(dataset_id, cluster)
            if ccc_data:
                outgoing = [e for e in ccc_data.get('edges', []) if e.get('from', '').endswith(cluster)]
                incoming = [e for e in ccc_data.get('edges', []) if e.get('to', '').endswith(cluster)]

                out_text = "\n".join([
                    f"- {e['ligand']} -> {e['receptor']} (to {e['to']}, pathway: {e.get('pathway', 'Unknown')})"
                    for e in outgoing[:10]
                ]) or "None detected"

                in_text = "\n".join([
                    f"- {e['ligand']} -> {e['receptor']} (from {e['from']}, pathway: {e.get('pathway', 'Unknown')})"
                    for e in incoming[:10]
                ]) or "None detected"

                sections.append(f"""## Cell-Cell Communication (Cluster {cluster})

### Outgoing Signals:
{out_text}

### Incoming Signals:
{in_text}""")
        except Exception:
            pass

    return "\n\n".join(sections)
