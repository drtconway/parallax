<script setup lang="ts">
import { computed, ref } from 'vue'
import { useAlignmentStore } from '../stores/alignment'

const store = useAlignmentStore()
const acceptError = ref<string | null>(null)

const result = computed(() => store.currentResult)
const readIndex = computed(() => store.currentResultIndex + 1)
const readTotal = computed(() => store.results.length)

const record = computed(() => {
  const r = result.value
  return r ? r.records[r.currentRecordIndex] : null
})
const alignmentIndex = computed(() => (result.value?.currentRecordIndex ?? 0) + 1)
const alignmentTotal = computed(() => result.value?.records.length ?? 0)

async function accept() {
  acceptError.value = null
  try {
    await store.acceptCurrent()
  } catch (e) {
    acceptError.value = e instanceof Error ? e.message : String(e)
  }
}
</script>

<template>
  <div v-if="result" class="record-navigator">
    <div class="read-nav">
      <button @click="store.prevResult" :disabled="store.currentResultIndex === 0">‹</button>
      <span>Read {{ readIndex }} / {{ readTotal }}</span>
      <button @click="store.nextResult" :disabled="store.currentResultIndex === readTotal - 1">›</button>
      <span class="status" :class="result.status">{{ result.status }}</span>
      <button
        v-if="result.status === 'pending'"
        @click="accept"
        class="accept-btn"
      >Accept</button>
    </div>

    <div class="alignment-nav" v-if="alignmentTotal > 0">
      <button @click="store.prevRecord" :disabled="alignmentIndex === 1">‹</button>
      <span>Alignment {{ alignmentIndex }} / {{ alignmentTotal }}</span>
      <button @click="store.nextRecord" :disabled="alignmentIndex === alignmentTotal">›</button>
      <span v-if="record" class="locus">{{ record.chrom }}:{{ record.start }}-{{ record.end }}</span>
    </div>

    <div v-if="acceptError" class="error">{{ acceptError }}</div>
  </div>
</template>
