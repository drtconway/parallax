import { defineStore } from 'pinia'
import { ref } from 'vue'

export interface AlignmentRecord {
  name: string
  chrom: string
  start: number
  end: number
}

export interface AlignmentResult {
  resultId: string
  bamUrl: string
  baiUrl: string
  contextTracks: ContextTrack[]
  records: AlignmentRecord[]
  currentIndex: number
}

export interface ContextTrack {
  name: string
  url: string
  indexURL: string
  format: string
}

export const useAlignmentStore = defineStore('alignment', () => {
  const result = ref<AlignmentResult | null>(null)
  const loading = ref(false)
  const error = ref<string | null>(null)

  async function align(fastqPath: string, contextBamPaths: string[]) {
    loading.value = true
    error.value = null
    result.value = null

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
      result.value = {
        resultId: data.result_id,
        bamUrl: data.bam_url,
        baiUrl: data.bai_url,
        contextTracks: data.context_tracks,
        records: data.records ?? [],
        currentIndex: 0,
      }
    } catch (e) {
      error.value = e instanceof Error ? e.message : String(e)
    } finally {
      loading.value = false
    }
  }

  function setRecords(records: AlignmentRecord[]) {
    if (result.value) {
      result.value.records = records
      result.value.currentIndex = 0
    }
  }

  function nextRecord() {
    if (result.value && result.value.currentIndex < result.value.records.length - 1) {
      result.value.currentIndex++
    }
  }

  function prevRecord() {
    if (result.value && result.value.currentIndex > 0) {
      result.value.currentIndex--
    }
  }

  return { result, loading, error, align, setRecords, nextRecord, prevRecord }
})
