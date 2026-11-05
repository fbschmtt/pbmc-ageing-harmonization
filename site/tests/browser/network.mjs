/** Keep browser journeys offline while preserving requests to their local test server. */
export async function installOfflineRoutes(page, appUrl) {
  const localOrigin = new URL(appUrl).origin
  await page.route('**/*', (route) => {
    const requestUrl = new URL(route.request().url())
    if (requestUrl.origin === localOrigin) return route.fallback()
    if (requestUrl.hostname === 'fonts.googleapis.com') {
      return route.fulfill({ status: 200, contentType: 'text/css', body: '' })
    }
    return route.abort()
  })
}
