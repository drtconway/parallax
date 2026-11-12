import { defineStore } from 'pinia'
import { ref, computed } from 'vue'

export interface AlignmentRecord {
  name: string
  chrom: string
  start: number
  end: number
}

export interface ContextTrack {
  name: string
  url: string
  indexURL: string
  format: string
}

export type ResultStatus = 'pending' | 'passing' | 'failing' | 'missing'

export interface ReadResult {
  digest: string
  status: ResultStatus
  bamUrl: string
  baiUrl: string
  expectedBamUrl: string | null
  expectedBaiUrl: string | null
  records: AlignmentRecord[]
  currentRecordIndex: number
}

export const useAlignmentStore = defineStore('alignment', () => {
  const results = ref<ReadResult[]>([])
  const currentResultIndex = ref(0)
  const contextTracks = ref<ContextTrack[]>([])
  const loading = ref(false)
  const error = ref<string | null>(null)

  const currentResult = computed(() => results.value[currentResultIndex.value] ?? null)

  // Only show pending and failing results — passing ones need no review
  const reviewableResults = computed(() =>
    results.value.filter(r => r.status === 'pending' || r.status === 'failing')
  )

  async function align(fastqPath: string, contextBamPaths: string[]) {
    loading.value = true
    error.value = null
    results.value = []
    currentResultIndex.value = 0

    try {
      const response = await fetch('/api/align', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ fastq_path: fastqPath, context_bam_paths: contextBamPaths }),
      })

      if (!response.ok) {
        const detail = await response.json().catch(() => ({ detail: response.statusText }))
        throw new Error(detail.detail ?? response.statusText)
      }

      const data = await response.json()
      contextTracks.value = data.context_tracks ?? []
      results.value = (data.results ?? []).map((r: any): ReadResult => ({
        digest: r.digest,
        status: r.status,
        bamUrl: r.bam_url,
        baiUrl: r.bai_url,
        expectedBamUrl: r.expected_bam_url ?? null,
        expectedBaiUrl: r.expected_bai_url ?? null,
        records: r.records ?? [],
        currentRecordIndex: 0,
      }))
      // Start at first reviewable result
      const firstReviewable = results.value.findIndex(r => r.status === 'pending' || r.status === 'failing')
      currentResultIndex.value = firstReviewable >= 0 ? firstReviewable : 0
    } catch (e) {
      error.value = e instanceof Error ? e.message : String(e)
    } finally {
      loading.value = false
    }
  }

  async function acceptCurrent() {
    const result = currentResult.value
    if (!result) return
    const response = await fetch(`/api/accept/${result.digest}`, { method: 'POST' })
    if (!response.ok) {
      const detail = await response.json().catch(() => ({ detail: response.statusText }))
      throw new Error(detail.detail ?? response.statusText)
    }
    result.status = 'passing'
    result.expectedBamUrl = result.bamUrl.replace('alignment.bam', 'expected.bam')
    result.expectedBaiUrl = result.baiUrl.replace('alignment.bam.bai', 'expected.bam.bai')
  }

  function nextResult() {
    if (currentResultIndex.value < results.value.length - 1) currentResultIndex.value++
  }

  function prevResult() {
    if (currentResultIndex.value > 0) currentResultIndex.value--
  }

  function nextRecord() {
    const r = currentResult.value
    if (r && r.currentRecordIndex < r.records.length - 1) r.currentRecordIndex++
  }

  function prevRecord() {
    const r = currentResult.value
    if (r && r.currentRecordIndex > 0) r.currentRecordIndex--
  }

  return {
    results, currentResultIndex, contextTracks, loading, error,
    currentResult, reviewableResults,
    align, acceptCurrent, nextResult, prevResult, nextRecord, prevRecord,
  }
})
