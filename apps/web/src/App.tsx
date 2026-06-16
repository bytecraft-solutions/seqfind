import { useState } from 'react'
import './App.css'

type Tool =
  | 'search'
  | 'stats'
  | 'transform'
  | 'subseq'
  | 'locate'
  | 'convert'

const SAMPLE_FASTA = `>chr1 demo record
ACGTATGCGTACGTTTATGCGTAATGCAT
>chr2 short example
TTTATGCGTAAATGCCGAT`

interface ToolMeta {
  id: Tool
  label: string
  description: string
}

const TOOLS: ToolMeta[] = [
  { id: 'search', label: 'search', description: 'Exact / wildcard / RC / approximate matching' },
  { id: 'stats', label: 'stats', description: 'N50, sum, GC%, Q20/Q30' },
  { id: 'transform', label: 'transform', description: 'case, RC, gap removal, DNA↔RNA, length filters' },
  { id: 'subseq', label: 'subseq', description: '1-based inclusive region extraction (negatives = from end)' },
  { id: 'locate', label: 'locate', description: 'IUPAC or regex motif localization' },
  { id: 'convert', label: 'convert', description: 'fq2fa, fx2tab, tab2fx' },
]

interface ToolParams {
  // search
  pattern: string
  reverse_complement: boolean
  max_edits: string
  // subseq
  region: string
  up_stream: string
  down_stream: string
  // locate
  patterns: string
  regex: boolean
  // transform
  uppercase: boolean
  lowercase: boolean
  remove_gaps: boolean
  reverse: boolean
  complement: boolean
  rc_full: boolean
  dna_to_rna: boolean
  rna_to_dna: boolean
  min_length: string
  max_length: string
  // convert
  convertKind: 'fq2fa' | 'fx2tab' | 'tab2fx'
}

const DEFAULTS: ToolParams = {
  pattern: 'ATGCGT',
  reverse_complement: false,
  max_edits: '',
  region: '1:10',
  up_stream: '0',
  down_stream: '0',
  patterns: 'ATG\nGCAT',
  regex: false,
  uppercase: false,
  lowercase: false,
  remove_gaps: false,
  reverse: false,
  complement: false,
  rc_full: false,
  dna_to_rna: false,
  rna_to_dna: false,
  min_length: '',
  max_length: '',
  convertKind: 'fx2tab',
}

async function runTool(tool: Tool, fasta: string, params: ToolParams): Promise<unknown> {
  const url = `/api/${tool}`
  let body: unknown
  switch (tool) {
    case 'search':
      body = {
        fasta,
        pattern: params.pattern,
        reverse_complement: params.reverse_complement,
        max_edits: params.max_edits ? Number(params.max_edits) : undefined,
      }
      break
    case 'stats':
      body = { fasta }
      break
    case 'transform':
      body = {
        fasta,
        uppercase: params.uppercase,
        lowercase: params.lowercase,
        remove_gaps: params.remove_gaps,
        reverse: params.reverse,
        complement: params.complement,
        reverse_complement: params.rc_full,
        dna_to_rna: params.dna_to_rna,
        rna_to_dna: params.rna_to_dna,
        min_length: params.min_length ? Number(params.min_length) : undefined,
        max_length: params.max_length ? Number(params.max_length) : undefined,
      }
      break
    case 'subseq':
      body = {
        fasta,
        region: params.region,
        up_stream: Number(params.up_stream) || 0,
        down_stream: Number(params.down_stream) || 0,
        reverse_complement: params.reverse_complement,
      }
      break
    case 'locate':
      body = {
        fasta,
        patterns: params.patterns.split('\n').map((s) => s.trim()).filter(Boolean),
        regex: params.regex,
        reverse_complement: params.reverse_complement,
      }
      break
    case 'convert':
      if (params.convertKind === 'tab2fx') {
        body = { kind: 'tab2fx', tab: fasta }
      } else {
        body = { kind: params.convertKind, fasta }
      }
      break
  }
  const resp = await fetch(url, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  })
  const text = await resp.text()
  if (!resp.ok) {
    throw new Error(`HTTP ${resp.status}: ${text}`)
  }
  try {
    return JSON.parse(text)
  } catch {
    return text
  }
}

function ParamsPanel({
  tool,
  params,
  setParams,
}: {
  tool: Tool
  params: ToolParams
  setParams: (p: ToolParams) => void
}) {
  const update = <K extends keyof ToolParams>(k: K, v: ToolParams[K]) =>
    setParams({ ...params, [k]: v })

  const text = (k: keyof ToolParams, label: string, placeholder = '') => (
    <label className="field">
      <span>{label}</span>
      <input
        type="text"
        value={String(params[k] as string)}
        placeholder={placeholder}
        onChange={(e) => update(k, e.target.value as never)}
      />
    </label>
  )

  const check = (k: keyof ToolParams, label: string) => (
    <label className="check">
      <input
        type="checkbox"
        checked={Boolean(params[k])}
        onChange={(e) => update(k, e.target.checked as never)}
      />
      <span>{label}</span>
    </label>
  )

  switch (tool) {
    case 'search':
      return (
        <div className="params">
          {text('pattern', 'Pattern', 'ATGCGT (? = wildcard)')}
          {text('max_edits', 'Max edits', 'e.g. 2 for Myers approx')}
          {check('reverse_complement', 'Search RC strand')}
        </div>
      )
    case 'stats':
      return <div className="params"><p className="hint">No parameters -POST FASTA and read the summary.</p></div>
    case 'transform':
      return (
        <div className="params">
          {check('uppercase', 'Uppercase')}
          {check('lowercase', 'Lowercase')}
          {check('remove_gaps', 'Remove gaps')}
          {check('reverse', 'Reverse')}
          {check('complement', 'Complement')}
          {check('rc_full', 'Reverse complement')}
          {check('dna_to_rna', 'DNA → RNA')}
          {check('rna_to_dna', 'RNA → DNA')}
          {text('min_length', 'Min length')}
          {text('max_length', 'Max length')}
        </div>
      )
    case 'subseq':
      return (
        <div className="params">
          {text('region', 'Region', 'START:END (1-based, negatives = from end)')}
          {text('up_stream', 'Up-stream context')}
          {text('down_stream', 'Down-stream context')}
          {check('reverse_complement', 'Take RC of slice')}
        </div>
      )
    case 'locate':
      return (
        <div className="params">
          <label className="field">
            <span>Patterns (one per line)</span>
            <textarea
              rows={4}
              value={params.patterns}
              onChange={(e) => update('patterns', e.target.value)}
            />
          </label>
          {check('regex', 'Treat as raw regex (else IUPAC)')}
          {check('reverse_complement', 'Also locate on RC strand')}
        </div>
      )
    case 'convert':
      return (
        <div className="params">
          <label className="field">
            <span>Direction</span>
            <select
              value={params.convertKind}
              onChange={(e) => update('convertKind', e.target.value as ToolParams['convertKind'])}
            >
              <option value="fq2fa">FASTQ → FASTA</option>
              <option value="fx2tab">Records → TSV (fx2tab)</option>
              <option value="tab2fx">TSV → records (tab2fx)</option>
            </select>
          </label>
          <p className="hint">
            For <code>tab2fx</code>, paste TSV in the input pane on the left.
          </p>
        </div>
      )
  }
}

function App() {
  const [tool, setTool] = useState<Tool>('search')
  const [fasta, setFasta] = useState(SAMPLE_FASTA)
  const [params, setParams] = useState<ToolParams>(DEFAULTS)
  const [result, setResult] = useState<string>('')
  const [error, setError] = useState<string>('')
  const [loading, setLoading] = useState(false)

  const run = async () => {
    setLoading(true)
    setError('')
    setResult('')
    try {
      const data = await runTool(tool, fasta, params)
      setResult(
        typeof data === 'string' ? data : JSON.stringify(data, null, 2),
      )
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e))
    } finally {
      setLoading(false)
    }
  }

  const currentTool = TOOLS.find((t) => t.id === tool)!

  return (
    <div className="app">
      <header className="topbar">
        <h1>seqfind workbench</h1>
        <a href="/api/docs" target="_blank" rel="noreferrer">
          API docs ↗
        </a>
      </header>

      <nav className="tabs">
        {TOOLS.map((t) => (
          <button
            key={t.id}
            className={t.id === tool ? 'tab active' : 'tab'}
            onClick={() => setTool(t.id)}
            type="button"
          >
            {t.label}
          </button>
        ))}
      </nav>

      <p className="description">{currentTool.description}</p>

      <main className="grid">
        <section className="pane">
          <h2>Input</h2>
          <textarea
            className="fasta"
            value={fasta}
            onChange={(e) => setFasta(e.target.value)}
            spellCheck={false}
          />
          <button type="button" onClick={() => setFasta(SAMPLE_FASTA)} className="ghost">
            Reset sample
          </button>
        </section>

        <section className="pane">
          <h2>Parameters</h2>
          <ParamsPanel tool={tool} params={params} setParams={setParams} />
          <button
            type="button"
            onClick={run}
            disabled={loading}
            className="primary"
          >
            {loading ? 'Running…' : `Run ${currentTool.label}`}
          </button>
        </section>

        <section className="pane full">
          <h2>Result</h2>
          {error && <pre className="error">{error}</pre>}
          <pre className="result">{result || '(no output yet)'}</pre>
        </section>
      </main>
    </div>
  )
}

export default App
