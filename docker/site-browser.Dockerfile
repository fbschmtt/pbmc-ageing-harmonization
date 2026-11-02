ARG SITE_IMAGE=pbmc-ageing-site:local
FROM ${SITE_IMAGE}

ENV PLAYWRIGHT_BROWSERS_PATH=/ms-playwright

RUN node node_modules/@playwright/test/cli.js install --with-deps chromium \
    && chmod -R a+rX /ms-playwright
