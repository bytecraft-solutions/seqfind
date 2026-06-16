//! PostgreSQL persistence for run history and results.
//!
//! Uses the blocking `postgres` crate (built on `tokio-postgres`) so it slots
//! into a synchronous CLI without an async runtime. The schema is created
//! lazily on first connection.

use anyhow::{Context, Result};
use postgres::{Client, NoTls};
use serde::Serialize;

use crate::matcher::Hit;
use crate::orf::Orf;
use crate::output::RecordStats;

pub const URL_ENV: &str = "SEQFIND_DB_URL";

/// Resolve the Postgres URL: explicit override > `SEQFIND_DB_URL` env > error.
pub fn resolve_url(override_url: Option<&str>) -> Result<String> {
    if let Some(url) = override_url.filter(|s| !s.is_empty()) {
        return Ok(url.to_string());
    }
    std::env::var(URL_ENV).map_err(|_| {
        anyhow::anyhow!(
            "{URL_ENV} is not set. Export it (e.g. \
             `export {URL_ENV}=postgresql://postgres:passw0rd@127.0.0.1/seqfind`) \
             or pass the URL explicitly with --save-db <URL> / --db <URL>."
        )
    })
}

const SCHEMA: &str = r#"
CREATE TABLE IF NOT EXISTS runs (
    id          BIGSERIAL PRIMARY KEY,
    ts          TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    subcommand  TEXT        NOT NULL,
    input_file  TEXT,
    params      JSONB       NOT NULL
);

CREATE TABLE IF NOT EXISTS hits (
    run_id   BIGINT NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    record   TEXT   NOT NULL,
    "start"  BIGINT NOT NULL,
    "end"    BIGINT NOT NULL,
    strand   TEXT   NOT NULL,
    matched  TEXT   NOT NULL,
    pattern  TEXT   NOT NULL
);

CREATE TABLE IF NOT EXISTS orfs (
    run_id      BIGINT NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    record      TEXT   NOT NULL,
    "start"     BIGINT NOT NULL,
    "end"       BIGINT NOT NULL,
    length      BIGINT NOT NULL,
    strand      TEXT   NOT NULL,
    start_codon TEXT   NOT NULL,
    stop_codon  TEXT   NOT NULL
);

CREATE TABLE IF NOT EXISTS record_stats (
    run_id     BIGINT NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    record     TEXT   NOT NULL,
    length     BIGINT NOT NULL,
    gc_percent DOUBLE PRECISION NOT NULL,
    a BIGINT NOT NULL,
    t BIGINT NOT NULL,
    g BIGINT NOT NULL,
    c BIGINT NOT NULL,
    n BIGINT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_hits_run  ON hits(run_id);
CREATE INDEX IF NOT EXISTS idx_orfs_run  ON orfs(run_id);
CREATE INDEX IF NOT EXISTS idx_stats_run ON record_stats(run_id);
"#;

pub struct Database {
    client: Client,
}

#[derive(Debug, Clone, Serialize)]
pub struct RunMeta {
    pub id: i64,
    pub ts: String,
    pub subcommand: String,
    pub input_file: Option<String>,
    pub params: String,
    pub result_count: i64,
}

impl Database {
    pub fn open(url: &str) -> Result<Self> {
        let mut client =
            Client::connect(url, NoTls).with_context(|| format!("failed to connect to {url}"))?;
        client
            .batch_execute(SCHEMA)
            .context("failed to apply schema")?;
        Ok(Self { client })
    }

    /// Insert a `runs` row and return its id.
    pub fn record_run(
        &mut self,
        subcommand: &str,
        input_file: Option<&str>,
        params_json: &str,
    ) -> Result<i64> {
        let params: serde_json::Value =
            serde_json::from_str(params_json).context("params is not valid JSON")?;
        let row = self.client.query_one(
            "INSERT INTO runs (subcommand, input_file, params)
             VALUES ($1, $2, $3) RETURNING id",
            &[&subcommand, &input_file, &params],
        )?;
        Ok(row.get(0))
    }

    pub fn insert_hits(&mut self, run_id: i64, hits: &[Hit]) -> Result<()> {
        if hits.is_empty() {
            return Ok(());
        }
        let mut tx = self.client.transaction()?;
        let stmt = tx.prepare(
            r#"INSERT INTO hits (run_id, record, "start", "end", strand, matched, pattern)
               VALUES ($1, $2, $3, $4, $5, $6, $7)"#,
        )?;
        for h in hits {
            tx.execute(
                &stmt,
                &[
                    &run_id,
                    &h.record,
                    &(h.start as i64),
                    &(h.end as i64),
                    &h.strand,
                    &h.matched,
                    &h.pattern,
                ],
            )?;
        }
        tx.commit()?;
        Ok(())
    }

    pub fn insert_orfs(&mut self, run_id: i64, orfs: &[Orf]) -> Result<()> {
        if orfs.is_empty() {
            return Ok(());
        }
        let mut tx = self.client.transaction()?;
        let stmt = tx.prepare(
            r#"INSERT INTO orfs (run_id, record, "start", "end", length, strand, start_codon, stop_codon)
               VALUES ($1, $2, $3, $4, $5, $6, $7, $8)"#,
        )?;
        for o in orfs {
            tx.execute(
                &stmt,
                &[
                    &run_id,
                    &o.record,
                    &(o.start as i64),
                    &(o.end as i64),
                    &(o.length as i64),
                    &o.strand,
                    &o.start_codon,
                    &o.stop_codon,
                ],
            )?;
        }
        tx.commit()?;
        Ok(())
    }

    pub fn insert_stats(&mut self, run_id: i64, stats: &[RecordStats]) -> Result<()> {
        if stats.is_empty() {
            return Ok(());
        }
        let mut tx = self.client.transaction()?;
        let stmt = tx.prepare(
            "INSERT INTO record_stats (run_id, record, length, gc_percent, a, t, g, c, n)
             VALUES ($1, $2, $3, $4, $5, $6, $7, $8, $9)",
        )?;
        for s in stats {
            tx.execute(
                &stmt,
                &[
                    &run_id,
                    &s.record,
                    &(s.length as i64),
                    &s.gc_percent,
                    &(s.a as i64),
                    &(s.t as i64),
                    &(s.g as i64),
                    &(s.c as i64),
                    &(s.n as i64),
                ],
            )?;
        }
        tx.commit()?;
        Ok(())
    }

    pub fn list_runs(&mut self, limit: usize) -> Result<Vec<RunMeta>> {
        let rows = self.client.query(
            "SELECT r.id,
                    to_char(r.ts AT TIME ZONE 'UTC', 'YYYY-MM-DD HH24:MI:SS'),
                    r.subcommand,
                    r.input_file,
                    r.params::text,
                    COALESCE((SELECT COUNT(*) FROM hits         WHERE run_id = r.id), 0) +
                    COALESCE((SELECT COUNT(*) FROM orfs         WHERE run_id = r.id), 0) +
                    COALESCE((SELECT COUNT(*) FROM record_stats WHERE run_id = r.id), 0) AS result_count
             FROM runs r
             ORDER BY r.id DESC
             LIMIT $1",
            &[&(limit as i64)],
        )?;
        Ok(rows
            .into_iter()
            .map(|row| RunMeta {
                id: row.get(0),
                ts: row.get(1),
                subcommand: row.get(2),
                input_file: row.get(3),
                params: row.get(4),
                result_count: row.get(5),
            })
            .collect())
    }

    pub fn get_run(&mut self, run_id: i64) -> Result<Option<RunMeta>> {
        let rows = self.client.query(
            "SELECT id,
                    to_char(ts AT TIME ZONE 'UTC', 'YYYY-MM-DD HH24:MI:SS'),
                    subcommand, input_file, params::text
             FROM runs WHERE id = $1",
            &[&run_id],
        )?;
        Ok(rows.into_iter().next().map(|row| RunMeta {
            id: row.get(0),
            ts: row.get(1),
            subcommand: row.get(2),
            input_file: row.get(3),
            params: row.get(4),
            result_count: 0,
        }))
    }

    pub fn get_hits(&mut self, run_id: i64) -> Result<Vec<Hit>> {
        let rows = self.client.query(
            r#"SELECT record, "start", "end", strand, matched, pattern
               FROM hits WHERE run_id = $1 ORDER BY record, "start""#,
            &[&run_id],
        )?;
        Ok(rows
            .into_iter()
            .map(|row| Hit {
                record: row.get(0),
                start: row.get::<_, i64>(1) as usize,
                end: row.get::<_, i64>(2) as usize,
                strand: row.get(3),
                matched: row.get(4),
                pattern: row.get(5),
            })
            .collect())
    }

    pub fn get_orfs(&mut self, run_id: i64) -> Result<Vec<Orf>> {
        let rows = self.client.query(
            r#"SELECT record, "start", "end", length, strand, start_codon, stop_codon
               FROM orfs WHERE run_id = $1 ORDER BY record, "start""#,
            &[&run_id],
        )?;
        Ok(rows
            .into_iter()
            .map(|row| Orf {
                record: row.get(0),
                start: row.get::<_, i64>(1) as usize,
                end: row.get::<_, i64>(2) as usize,
                length: row.get::<_, i64>(3) as usize,
                strand: row.get(4),
                start_codon: row.get(5),
                stop_codon: row.get(6),
            })
            .collect())
    }

    pub fn get_stats(&mut self, run_id: i64) -> Result<Vec<RecordStats>> {
        let rows = self.client.query(
            "SELECT record, length, gc_percent, a, t, g, c, n
             FROM record_stats WHERE run_id = $1 ORDER BY record",
            &[&run_id],
        )?;
        Ok(rows
            .into_iter()
            .map(|row| RecordStats {
                record: row.get(0),
                length: row.get::<_, i64>(1) as usize,
                gc_percent: row.get(2),
                a: row.get::<_, i64>(3) as usize,
                t: row.get::<_, i64>(4) as usize,
                g: row.get::<_, i64>(5) as usize,
                c: row.get::<_, i64>(6) as usize,
                n: row.get::<_, i64>(7) as usize,
            })
            .collect())
    }
}
