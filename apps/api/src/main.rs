use axum::{Json, http::StatusCode, response::IntoResponse};
use serde::{Deserialize, Serialize};
use std::net::SocketAddr;
use utoipa::{OpenApi, ToSchema};
use utoipa_axum::{router::OpenApiRouter, routes};
use utoipa_scalar::{Scalar, Servable};

use seqfind_core::{matcher, seq};

#[derive(OpenApi)]
#[openapi(
    info(
        title = "seqfind-api",
        version = "0.1.0",
        description = "HTTP API over the seqfind-core deterministic search library.",
    ),
    tags(
        (name = "search", description = "Pattern search over FASTA"),
        (name = "meta",   description = "Service metadata"),
    ),
)]
struct ApiDoc;

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt::init();

    let (router, api) = OpenApiRouter::with_openapi(ApiDoc::openapi())
        .routes(routes!(health))
        .routes(routes!(search))
        .split_for_parts();

    let app = router.merge(Scalar::with_url("/docs", api.clone())).route(
        "/openapi.json",
        axum::routing::get(move || {
            let api = api.clone();
            async move { Json(api) }
        }),
    );

    let addr: SocketAddr = "0.0.0.0:8080".parse()?;
    tracing::info!("seqfind-api listening on http://{addr}  (docs: /docs)");
    let listener = tokio::net::TcpListener::bind(addr).await?;
    axum::serve(listener, app).await?;
    Ok(())
}

#[utoipa::path(get, path = "/health", tag = "meta", responses((status = 200, body = String)))]
async fn health() -> &'static str {
    "ok"
}

#[derive(Debug, Deserialize, ToSchema)]
pub struct SearchRequest {
    /// FASTA payload (one or more records, plain text).
    #[schema(example = ">chr1\nACGTATGCGTACG\n>chr2\nTTTATGCGTAA")]
    pub fasta: String,
    /// Pattern to search for. `?` is a single-base wildcard.
    #[schema(example = "ATGCGT")]
    pub pattern: String,
    /// Also search the reverse complement and report `-` strand hits.
    #[serde(default)]
    pub reverse_complement: bool,
    /// Allow up to N edits (Myers approximate matching). 0 = exact.
    #[serde(default)]
    pub max_edits: Option<usize>,
}

#[derive(Debug, Serialize, ToSchema)]
pub struct SearchHit {
    pub record: String,
    pub start: usize,
    pub end: usize,
    pub strand: String,
    pub matched: String,
    pub pattern: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub edits: Option<usize>,
}

impl From<matcher::Hit> for SearchHit {
    fn from(h: matcher::Hit) -> Self {
        Self {
            record: h.record,
            start: h.start,
            end: h.end,
            strand: h.strand,
            matched: h.matched,
            pattern: h.pattern,
            edits: h.edits,
        }
    }
}

#[derive(Debug, Serialize, ToSchema)]
pub struct SearchResponse {
    pub count: usize,
    pub hits: Vec<SearchHit>,
}

#[derive(Debug, Serialize, ToSchema)]
pub struct ApiError {
    pub error: String,
}

/// Search a FASTA payload for a pattern (with optional wildcard, RC, or
/// approximate matching). Returns hits in the same shape as the CLI.
#[utoipa::path(
    post,
    path = "/search",
    tag = "search",
    request_body = SearchRequest,
    responses(
        (status = 200, body = SearchResponse),
        (status = 400, body = ApiError, description = "Invalid FASTA or empty pattern"),
    ),
)]
async fn search(Json(req): Json<SearchRequest>) -> Result<Json<SearchResponse>, ApiErrorResponse> {
    if req.pattern.is_empty() {
        return Err(ApiErrorResponse::bad_request("pattern must not be empty"));
    }

    let records =
        seq::read_seqs_str(&req.fasta).map_err(|e| ApiErrorResponse::bad_request(e.to_string()))?;

    let mut hits = Vec::new();
    for r in &records {
        let found = match req.max_edits {
            Some(edits) if edits > 0 => matcher::find_matches_with_max_edits(
                &r.name,
                &r.sequence,
                &req.pattern,
                req.reverse_complement,
                edits,
            ),
            _ => matcher::find_matches(&r.name, &r.sequence, &req.pattern, req.reverse_complement),
        };
        hits.extend(found.into_iter().map(SearchHit::from));
    }

    Ok(Json(SearchResponse {
        count: hits.len(),
        hits,
    }))
}

struct ApiErrorResponse {
    status: StatusCode,
    body: ApiError,
}

impl ApiErrorResponse {
    fn bad_request(msg: impl Into<String>) -> Self {
        Self {
            status: StatusCode::BAD_REQUEST,
            body: ApiError { error: msg.into() },
        }
    }
}

impl IntoResponse for ApiErrorResponse {
    fn into_response(self) -> axum::response::Response {
        (self.status, Json(self.body)).into_response()
    }
}
