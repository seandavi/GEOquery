"readUrl" <-
function (url) 
{
    options(show.error.messages = FALSE)
    con <- try(url(url, open = "r"))
    options(show.error.messages = TRUE)
    if (inherits(con, "try-error")) {
        stop(paste("Can't connent to url", url))
    }
    temp <- readLines(con)
    close(con)
    return(temp)
}
