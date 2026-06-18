<script setup lang="ts">
import { computed } from 'vue'
import { useAlignmentStore } from '../stores/alignment'

const store = useAlignmentStore()

const current = computed(() => store.result?.records[store.result.currentIndex])
const total = computed(() => store.result?.records.length ?? 0)
const index = computed(() => (store.result?.currentIndex ?? 0) + 1)
</script>

<template>
  <div v-if="total > 0" class="record-navigator">
    <button @click="store.prevRecord" :disabled="index === 1">‹</button>
    <span>{{ index }} / {{ total }}</span>
    <button @click="store.nextRecord" :disabled="index === total">›</button>
    <span v-if="current" class="locus">{{ current.chrom }}:{{ current.start }}-{{ current.end }}</span>
  </div>
</template>
