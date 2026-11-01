// Parse RFC-style quoted fields, including embedded commas and newlines.
export function parseCsv(text) {
  const rows = []
  let row = [],
    field = '',
    quoted = false,
    closed = false,
    atStart = true
  const pushField = () => {
    row.push(field)
    field = ''
    closed = false
    atStart = true
  }
  const pushRow = () => {
    pushField()
    rows.push(row)
    row = []
  }
  for (let i = 0; i < text.length; i++) {
    const ch = text[i]
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') {
        field += '"'
        i++
      } else if (ch === '"') {
        quoted = false
        closed = true
      } else field += ch
    } else if (ch === '"') {
      if (!atStart || closed) throw new Error('Invalid quote in CSV field')
      quoted = true
      atStart = false
    } else if (ch === ',') pushField()
    else if (ch === '\n' || ch === '\r') {
      if (ch === '\r' && text[i + 1] === '\n') i++
      pushRow()
    } else {
      if (closed) throw new Error('Unexpected content after quoted CSV field')
      field += ch
      atStart = false
    }
  }
  if (quoted) throw new Error('Unterminated quoted CSV field')
  if (row.length || field.length || closed) pushRow()
  const headers = rows.shift() || []
  if (!headers.length || headers.some((h) => !h) || new Set(headers).size !== headers.length)
    throw new Error('Invalid CSV headers')
  if (rows.some((r) => r.length !== headers.length)) throw new Error('CSV row has wrong number of fields')
  return { headers, rows }
}
