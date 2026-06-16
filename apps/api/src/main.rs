use axum::{Router, routing::get};
use std::net::SocketAddr;

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt::init();

    let app = Router::new().route("/health", get(health));

    let addr: SocketAddr = "0.0.0.0:8080".parse()?;
    tracing::info!("seqfind-api listening on http://{addr}");
    let listener = tokio::net::TcpListener::bind(addr).await?;
    axum::serve(listener, app).await?;
    Ok(())
}

async fn health() -> &'static str {
    "ok"
}
