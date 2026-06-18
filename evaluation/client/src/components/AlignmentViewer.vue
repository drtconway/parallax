<script setup lang="ts">
import { ref, watch, onMounted, onUnmounted, nextTick } from 'vue'
import igv from 'igv'
import type { Browser } from 'igv'
import { useAlignmentStore } from '../stores/alignment'
import type { AlignmentRecord } from '../stores/alignment'

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

async function initBrowser(bamUrl: string, baiUrl: string, contextTracks: any[]) {
  if (!container.value) return

  if (browser) {
    igv.removeBrowser(browser)
    browser = null
  }

  const tracks = [
    {
      name: 'Alignment',
      url: bamUrl,
      indexURL: baiUrl,
      format: 'bam',
      type: 'alignment',
    },
    ...contextTracks.map((t: any) => ({ ...t, type: 'alignment' })),
  ]

  browser = await igv.createBrowser(container.value, {
    genome: 'hg38',
    tracks,
  })

  const records = store.result?.records
  if (records && records.length > 0) {
    browser.search(expandedLocus(records[0]))
  }
}

watch(
  () => store.result?.currentIndex,
  (idx) => {
    if (!browser || !store.result || idx == null) return
    const rec = store.result.records[idx]
    if (rec) browser.search(expandedLocus(rec))
  }
)

onMounted(() => {
  if (store.result) {
    nextTick(() => initBrowser(store.result!.bamUrl, store.result!.baiUrl, store.result!.contextTracks))
  }
})

watch(
  () => store.result?.resultId,
  () => {
    if (store.result) {
      nextTick(() => initBrowser(store.result!.bamUrl, store.result!.baiUrl, store.result!.contextTracks))
    }
  }
)

onUnmounted(() => {
  if (browser) igv.removeBrowser(browser)
})
</script>

<template>
  <div ref="container" class="igv-container"></div>
</template>
