<script setup lang="ts">
import { ref, watch, onMounted, onUnmounted, nextTick } from 'vue'
import igv from 'igv'
import type { Browser } from 'igv'
import { useAlignmentStore } from '../stores/alignment'
import type { AlignmentRecord, ReadResult } from '../stores/alignment'

const store = useAlignmentStore()
const container = ref<HTMLElement | null>(null)
let browser: Browser | null = null

const LOCUS_PADDING_DIVISOR = 6
const LOCUS_MIN_WIDTH = 50

function expandedLocus(rec: AlignmentRecord): string {
  const len = rec.end - rec.start
  const padding = Math.max(Math.round(len / LOCUS_PADDING_DIVISOR), Math.round((LOCUS_MIN_WIDTH - len) / 2 + 1))
  return `${rec.chrom}:${Math.max(1, rec.start - padding)}-${rec.end + padding}`
}

async function initBrowser(result: ReadResult) {
  if (!container.value) return

  if (browser) {
    igv.removeBrowser(browser)
    browser = null
  }

  const tracks: any[] = [
    {
      name: 'Alignment',
      url: result.bamUrl,
      indexURL: result.baiUrl,
      format: 'bam',
      type: 'alignment',
    },
  ]

  if (result.expectedBamUrl && result.expectedBaiUrl) {
    tracks.push({
      name: 'Expected',
      url: result.expectedBamUrl,
      indexURL: result.expectedBaiUrl,
      format: 'bam',
      type: 'alignment',
    })
  }

  for (const t of store.contextTracks) {
    tracks.push({ ...t, type: 'alignment' })
  }

  browser = await igv.createBrowser(container.value, {
    genome: 'hg38',
    tracks,
  })

  if (result.records.length > 0) {
    browser.search(expandedLocus(result.records[0]))
  }
}

watch(
  () => store.currentResult?.currentRecordIndex,
  (idx) => {
    if (!browser || !store.currentResult || idx == null) return
    const rec = store.currentResult.records[idx]
    if (rec) browser.search(expandedLocus(rec))
  }
)

watch(
  () => store.currentResult?.digest,
  () => {
    if (store.currentResult) {
      nextTick(() => initBrowser(store.currentResult!))
    }
  }
)

onMounted(() => {
  if (store.currentResult) {
    nextTick(() => initBrowser(store.currentResult!))
  }
})

onUnmounted(() => {
  if (browser) igv.removeBrowser(browser)
})
</script>

<template>
  <div ref="container" class="igv-container"></div>
</template>

<style scoped>
.igv-container {
  width: 100%;
  font-size: 12px;
  line-height: normal;
}
</style>
