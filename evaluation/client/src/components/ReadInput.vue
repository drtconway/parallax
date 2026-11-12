<script setup lang="ts">
import { ref } from 'vue'
import { useAlignmentStore } from '../stores/alignment'

const store = useAlignmentStore()
const fastqPath = ref('')
const contextPaths = ref<string[]>([''])

function addContextPath() {
  contextPaths.value.push('')
}

function removeContextPath(i: number) {
  contextPaths.value.splice(i, 1)
}

async function submit() {
  const paths = contextPaths.value.filter(p => p.trim().length > 0)
  await store.align(fastqPath.value.trim(), paths)
}
</script>

<template>
  <div class="read-input">
    <div class="field">
      <label>FASTQ path</label>
      <input v-model="fastqPath" placeholder="/path/to/reads.fastq.gz" />
    </div>

    <div class="field">
      <label>Context BAM files</label>
      <div v-for="(_, i) in contextPaths" :key="i" class="context-row">
        <input v-model="contextPaths[i]" placeholder="/path/to/context.bam" />
        <button @click="removeContextPath(i)" :disabled="contextPaths.length === 1">✕</button>
      </div>
      <button @click="addContextPath">+ Add BAM</button>
    </div>

    <button @click="submit" :disabled="!fastqPath || store.loading">
      {{ store.loading ? 'Aligning…' : 'Align' }}
    </button>

    <div v-if="store.error" class="error">{{ store.error }}</div>
  </div>
</template>
