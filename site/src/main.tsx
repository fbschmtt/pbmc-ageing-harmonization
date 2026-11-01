import { createRoot } from 'react-dom/client'
import { App } from './app/App'
import { installVisualEncoding } from './theme/visualEncoding'
import './styles.css'

installVisualEncoding(document.documentElement)
createRoot(document.getElementById('root')!).render(<App />)
