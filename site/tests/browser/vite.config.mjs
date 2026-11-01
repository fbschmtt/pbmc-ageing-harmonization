import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({ plugins: [react()], publicDir: process.env.ATLAS_TEST_PUBLIC_DIR })
