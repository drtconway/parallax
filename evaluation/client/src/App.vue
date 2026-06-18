<script setup lang="ts">
import ReadInput from './components/ReadInput.vue'
import AlignmentViewer from './components/AlignmentViewer.vue'
import RecordNavigator from './components/RecordNavigator.vue'
import { useAlignmentStore } from './stores/alignment'

const store = useAlignmentStore()
</script>

<template>
  <div class="app">
    <header>
      <h1>Parallax Evaluation</h1>
    </header>
    <main>
      <ReadInput />
      <RecordNavigator />
      <AlignmentViewer v-if="store.currentResult" />
    </main>
  </div>
</template>
