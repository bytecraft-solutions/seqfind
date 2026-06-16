use axum::{Json, http::StatusCode, response::IntoResponse};
use serde::{Deserialize, Serialize};
use std::net::SocketAddr;
use utoipa::{OpenApi, ToSchema};
use utoipa_axum::{router::OpenApiRouter, routes};
use utoipa_scalar::{Scalar, Servable};

use seqfind_core::{
    convert as core_convert, locate, matcher, seq, subseq as core_subseq, summary, transform,
};

#[derive(OpenApi)]
#[openapi(
    info(
        title = "seqfind-api",
        version = "0.1.0",
        description = "HTTP API over the seqfind-core deterministic search library.",
    ),
    tags(
        (name = "search",    description = "Pattern search over FASTA"),
        (name = "stats",     description = "Aggregate sequence statistics"),
        (name = "transform", description = "Per-record sequence transforms"),
        (name = "subseq",    description = "Region extraction"),
        (name = "locate",    description = "Motif / IUPAC / regex localization"),
        (name = "convert",   description = "Format converters (fq2fa, fx2tab, tab2fx)"),
        (name = "meta",      description = "Service metadata"),
    ),
)]
struct ApiDoc;

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt::init();

    let (router, api) = OpenApiRouter::with_openapi(ApiDoc::openapi())
        .routes(routes!(health))
        .routes(routes!(search))
        .routes(routes!(stats))
        .routes(routes!(transform_endpoint))
        .routes(routes!(subseq_endpoint))
        .routes(routes!(locate_endpoint))
        .routes(routes!(convert_endpoint))
        .split_for_parts();

    let app = router
        .merge(Scalar::with_url("/docs", api.clone()))
        .route(
            "/openapi.json",
            axum::routing::get(move || {
                let api = api.clone();
                async move { Json(api) }
            }),
        )
        .layer(tower_http::cors::CorsLayer::permissive());

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

// ---------- /search ----------

#[derive(Debug, Deserialize, ToSchema)]
pub struct SearchRequest {
    #[schema(example = ">chr1\nACGTATGCGTACG\n>chr2\nTTTATGCGTAA")]
    pub fasta: String,
    #[schema(example = "ATGCGT")]
    pub pattern: String,
    #[serde(default)]
    pub reverse_complement: bool,
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

#[utoipa::path(
    post, path = "/search", tag = "search",
    request_body = SearchRequest,
    responses(
        (status = 200, body = SearchResponse),
        (status = 400, body = ApiError),
    ),
)]
async fn search(Json(req): Json<SearchRequest>) -> Result<Json<SearchResponse>, ApiErrorResponse> {
    if req.pattern.is_empty() {
        return Err(ApiErrorResponse::bad_request("pattern must not be empty"));
    }
    let records = parse_fasta(&req.fasta)?;
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

// ---------- /stats ----------

#[derive(Debug, Deserialize, ToSchema)]
pub struct StatsRequest {
    #[schema(example = ">chr1\nACGTATGCGTACG")]
    pub fasta: String,
}

#[derive(Debug, Serialize, ToSchema)]
pub struct StatsResponse {
    pub num_seqs: usize,
    pub sum_len: usize,
    pub min_len: usize,
    pub avg_len: f64,
    pub max_len: usize,
    pub n50: usize,
    pub gc_percent: f64,
    pub n_count: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub q20_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub q30_percent: Option<f64>,
}

impl From<summary::StatsSummary> for StatsResponse {
    fn from(s: summary::StatsSummary) -> Self {
        Self {
            num_seqs: s.num_seqs,
            sum_len: s.sum_len,
            min_len: s.min_len,
            avg_len: s.avg_len,
            max_len: s.max_len,
            n50: s.n50,
            gc_percent: s.gc_percent,
            n_count: s.n_count,
            q20_percent: s.q20_percent,
            q30_percent: s.q30_percent,
        }
    }
}

#[utoipa::path(
    post, path = "/stats", tag = "stats",
    request_body = StatsRequest,
    responses((status = 200, body = StatsResponse), (status = 400, body = ApiError)),
)]
async fn stats(Json(req): Json<StatsRequest>) -> Result<Json<StatsResponse>, ApiErrorResponse> {
    let records = parse_fasta(&req.fasta)?;
    Ok(Json(summary::compute(&records).into()))
}

// ---------- /transform ----------

#[derive(Debug, Deserialize, ToSchema)]
pub struct TransformRequest {
    pub fasta: String,
    #[serde(default)]
    pub uppercase: bool,
    #[serde(default)]
    pub lowercase: bool,
    #[serde(default)]
    pub remove_gaps: bool,
    #[serde(default)]
    pub dna_to_rna: bool,
    #[serde(default)]
    pub rna_to_dna: bool,
    #[serde(default)]
    pub reverse: bool,
    #[serde(default)]
    pub complement: bool,
    #[serde(default)]
    pub reverse_complement: bool,
    #[serde(default)]
    pub min_length: Option<usize>,
    #[serde(default)]
    pub max_length: Option<usize>,
}

#[derive(Debug, Serialize, ToSchema)]
pub struct RecordOut {
    pub name: String,
    pub sequence: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub quality: Option<String>,
}

impl From<seq::SeqRecord> for RecordOut {
    fn from(r: seq::SeqRecord) -> Self {
        Self {
            name: r.name,
            sequence: r.sequence,
            quality: r.quality.map(|q| String::from_utf8_lossy(&q).into_owned()),
        }
    }
}

#[derive(Debug, Serialize, ToSchema)]
pub struct RecordsResponse {
    pub count: usize,
    pub records: Vec<RecordOut>,
}

#[utoipa::path(
    post, path = "/transform", tag = "transform",
    request_body = TransformRequest,
    responses((status = 200, body = RecordsResponse), (status = 400, body = ApiError)),
)]
async fn transform_endpoint(
    Json(req): Json<TransformRequest>,
) -> Result<Json<RecordsResponse>, ApiErrorResponse> {
    let records = parse_fasta(&req.fasta)?;
    let opts = transform::TransformOptions {
        uppercase: req.uppercase,
        lowercase: req.lowercase,
        remove_gaps: req.remove_gaps,
        dna_to_rna: req.dna_to_rna,
        rna_to_dna: req.rna_to_dna,
        reverse: req.reverse,
        complement: req.complement,
        reverse_complement: req.reverse_complement,
        min_length: req.min_length,
        max_length: req.max_length,
    };
    let out: Vec<RecordOut> = records
        .iter()
        .filter_map(|r| transform::apply(r, &opts))
        .map(RecordOut::from)
        .collect();
    Ok(Json(RecordsResponse {
        count: out.len(),
        records: out,
    }))
}

// ---------- /subseq ----------

#[derive(Debug, Deserialize, ToSchema)]
pub struct SubseqRequest {
    pub fasta: String,
    /// `START:END` (1-based inclusive, negatives count from end).
    #[schema(example = "1:10")]
    pub region: String,
    #[serde(default)]
    pub up_stream: usize,
    #[serde(default)]
    pub down_stream: usize,
    #[serde(default)]
    pub reverse_complement: bool,
}

#[utoipa::path(
    post, path = "/subseq", tag = "subseq",
    request_body = SubseqRequest,
    responses((status = 200, body = RecordsResponse), (status = 400, body = ApiError)),
)]
async fn subseq_endpoint(
    Json(req): Json<SubseqRequest>,
) -> Result<Json<RecordsResponse>, ApiErrorResponse> {
    let region = core_subseq::Region::parse(&req.region)
        .map_err(|e| ApiErrorResponse::bad_request(e.to_string()))?;
    let opts = core_subseq::SubseqOptions {
        region,
        reverse_complement: req.reverse_complement,
        up_stream: req.up_stream,
        down_stream: req.down_stream,
    };
    let records = parse_fasta(&req.fasta)?;
    let mut out = Vec::new();
    for r in &records {
        let extracted = core_subseq::subseq_record(r, &opts)
            .map_err(|e| ApiErrorResponse::bad_request(e.to_string()))?;
        if let Some(rec) = extracted {
            out.push(RecordOut::from(rec));
        }
    }
    Ok(Json(RecordsResponse {
        count: out.len(),
        records: out,
    }))
}

// ---------- /locate ----------

#[derive(Debug, Deserialize, ToSchema)]
pub struct LocateRequest {
    pub fasta: String,
    pub patterns: Vec<String>,
    #[serde(default)]
    pub regex: bool,
    #[serde(default)]
    pub reverse_complement: bool,
}

#[derive(Debug, Serialize, ToSchema)]
pub struct LocateHitOut {
    pub record: String,
    pub pattern_name: String,
    pub pattern: String,
    pub start: usize,
    pub end: usize,
    pub strand: String,
    pub matched: String,
}

impl From<locate::LocateHit> for LocateHitOut {
    fn from(h: locate::LocateHit) -> Self {
        Self {
            record: h.record,
            pattern_name: h.pattern_name,
            pattern: h.pattern,
            start: h.start,
            end: h.end,
            strand: h.strand,
            matched: h.matched,
        }
    }
}

#[derive(Debug, Serialize, ToSchema)]
pub struct LocateResponse {
    pub count: usize,
    pub hits: Vec<LocateHitOut>,
}

#[utoipa::path(
    post, path = "/locate", tag = "locate",
    request_body = LocateRequest,
    responses((status = 200, body = LocateResponse), (status = 400, body = ApiError)),
)]
async fn locate_endpoint(
    Json(req): Json<LocateRequest>,
) -> Result<Json<LocateResponse>, ApiErrorResponse> {
    if req.patterns.is_empty() {
        return Err(ApiErrorResponse::bad_request(
            "at least one pattern required",
        ));
    }
    let patterns: Vec<locate::LocatePattern> = req
        .patterns
        .iter()
        .enumerate()
        .map(|(i, raw)| {
            let name = format!("p{}", i + 1);
            if req.regex {
                locate::LocatePattern::new(name, raw.clone())
            } else {
                locate::LocatePattern::from_iupac(name, raw.clone())
            }
        })
        .collect::<Result<_, _>>()
        .map_err(|e| ApiErrorResponse::bad_request(e.to_string()))?;
    let records = parse_fasta(&req.fasta)?;
    let mut hits = Vec::new();
    for r in &records {
        hits.extend(
            locate::locate_in_record(r, &patterns, req.reverse_complement)
                .into_iter()
                .map(LocateHitOut::from),
        );
    }
    Ok(Json(LocateResponse {
        count: hits.len(),
        hits,
    }))
}

// ---------- /convert ----------

#[derive(Debug, Deserialize, ToSchema)]
#[serde(rename_all = "lowercase", tag = "kind")]
pub enum ConvertRequest {
    /// FASTQ → FASTA (drop quality). Input body: FASTQ in `fasta`.
    Fq2fa { fasta: String },
    /// Records → TSV. Input body: FASTA or FASTQ in `fasta`.
    Fx2tab { fasta: String },
    /// TSV → records. Input body: TSV in `tab`.
    Tab2fx { tab: String },
}

#[derive(Debug, Serialize, ToSchema)]
pub struct ConvertResponse {
    /// Plain-text payload in the target format (FASTA/TSV/records).
    pub output: String,
}

#[utoipa::path(
    post, path = "/convert", tag = "convert",
    request_body = ConvertRequest,
    responses((status = 200, body = ConvertResponse), (status = 400, body = ApiError)),
)]
async fn convert_endpoint(
    Json(req): Json<ConvertRequest>,
) -> Result<Json<ConvertResponse>, ApiErrorResponse> {
    match req {
        ConvertRequest::Fq2fa { fasta } => {
            let records = parse_fasta(&fasta)?;
            let converted = core_convert::fq2fa(&records);
            let mut buf = Vec::new();
            for r in &converted {
                seq::write_fasta_record(&mut buf, &r.name, &r.sequence, 70)
                    .map_err(|e| ApiErrorResponse::bad_request(e.to_string()))?;
            }
            Ok(Json(ConvertResponse {
                output: String::from_utf8(buf).unwrap_or_default(),
            }))
        }
        ConvertRequest::Fx2tab { fasta } => {
            let records = parse_fasta(&fasta)?;
            let mut buf = Vec::new();
            core_convert::fx2tab(&mut buf, &records)
                .map_err(|e| ApiErrorResponse::bad_request(e.to_string()))?;
            Ok(Json(ConvertResponse {
                output: String::from_utf8(buf).unwrap_or_default(),
            }))
        }
        ConvertRequest::Tab2fx { tab } => {
            let records = core_convert::tab2fx(std::io::Cursor::new(tab.as_bytes().to_vec()))
                .map_err(|e| ApiErrorResponse::bad_request(e.to_string()))?;
            let mut buf = Vec::new();
            let any_qual = records.iter().any(|r| r.quality.is_some());
            for r in &records {
                if any_qual {
                    seq::write_fastq_record(&mut buf, r)
                        .map_err(|e| ApiErrorResponse::bad_request(e.to_string()))?;
                } else {
                    seq::write_fasta_record(&mut buf, &r.name, &r.sequence, 70)
                        .map_err(|e| ApiErrorResponse::bad_request(e.to_string()))?;
                }
            }
            Ok(Json(ConvertResponse {
                output: String::from_utf8(buf).unwrap_or_default(),
            }))
        }
    }
}

// ---------- shared error helpers ----------

#[derive(Debug, Serialize, ToSchema)]
pub struct ApiError {
    pub error: String,
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

fn parse_fasta(s: &str) -> Result<Vec<seq::SeqRecord>, ApiErrorResponse> {
    seq::read_seqs_str(s).map_err(|e| ApiErrorResponse::bad_request(e.to_string()))
}
