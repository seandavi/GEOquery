# ADR-0007: Access private GEO records via a reviewer token

- Status: accepted
- Date: 2026-07-17
- Deciders: Sean Davis

## Context

NCBI GEO lets authors share a private/embargoed Series with reviewers through a
**reviewer access token** — an anonymous, read-only, expiring credential created
from the "Reviewer access" link on the GSE record. Users with such a token had
no way to load the data with GEOquery (issue #154); `getGEO()` returned HTTP 404
for private accessions.

Two facts constrain the design:

1. **Private records are not on the GEO FTP tree.** The FTP paths GEOquery uses
   for the fast Series Matrix path (`ftp.ncbi.nlm.nih.gov/geo/series/<stub>/<GSE>/matrix/`)
   and for the family SOFT (`.../soft/<GSE>_family.soft.gz`) return 404 for a
   private Series. Only the CGI endpoint `acc.cgi` honors the token.
2. **`acc.cgi` returns SOFT text, not a Series Matrix.** So a token-authenticated
   GSE fetch must go through the SOFT parser (`parseGSE`) and yields a `GSE` S4
   object, not an `ExpressionSet`/`SummarizedExperiment`.

The token is supplied to `acc.cgi` as a `&token=<TOKEN>` query parameter, e.g.
`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123456&targ=all&form=text&view=full&token=<TOKEN>`.

## Decision

We will add an optional `token` argument to `getGEO()` and `getGEOfile()`.

- When `token` is supplied and the accession is a GSE, `getGEO()` bypasses the
  FTP Series Matrix path (`getAndParseGSEMatrices()`) and uses the `getGEOfile()`
  SOFT path, returning a `GSE` S4 object. This is documented as a deliberate
  return-type difference (the same shape as `GSEMatrix = FALSE`).
- `getGEOfile()` appends `token=<token>` to the `acc.cgi` request for GSE
  (brief/quick and, with a token, the full family via `targ=all`), GPL, and GSM.
  FTP URLs (GDS, Annotation-GPL) never carry a token.
- A single shared internal helper, `.append_token()`, adds the query parameter,
  keeping the URL construction DRY.

## Consequences

- Users with a reviewer token can load private Series (as `GSE` objects) and
  private GSM/GPL records without leaving GEOquery.
- The token lives only in the request URL; no change to `downloadFile()` /
  `.geo_request()` is needed, and the existing httr2 mockability (`.geo_request()`
  returns an unperformed request; ADR/#173) lets us unit-test URL construction.
- We cannot exercise a real token in CI, so automated coverage is limited to
  URL-construction unit tests (token appended; `acc.cgi` path forced for a
  token'd GSE). The end-to-end path is validated by user report.
- Tokens are secrets: they appear in URLs and therefore potentially in logs and
  the on-disk cache key. We document that tokens are sensitive and short-lived;
  we do not persist them beyond the call.

## Alternatives considered

- **Support a token only on the FTP path.** Rejected: private records are simply
  absent from FTP, so there is nothing to authenticate against there.
- **Add a general HTTP header / auth mechanism.** Rejected as over-engineered;
  GEO's reviewer access is specifically a URL query parameter, and a bespoke
  `token` argument is the least surprising interface.
- **Return a Series Matrix for token'd GSEs by scraping.** Rejected: `acc.cgi`
  yields SOFT, and reusing the existing SOFT parser is simpler and correct; the
  return-type difference is documented.
