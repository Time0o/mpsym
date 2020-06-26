transformColumns <- function(d, column, f, ...) {
    for (i in 1:length(d[[column]])) {
        s <- d[[column]][i]

        s <- f(s)
        for (f_ in list(...))
            s <- f_(s)

        d[[column]][i] <- s
    }

    d
}

stripPrefix <- function(s) {
    gsub('^.*? ', '', s)
}

translate <- function(dict) {
    function(s) {
        if (!is.null(dict[[s]]))
            s <- dict[[s]]

        s
    }
}

capitalize <- function(s) {
    s1 <- substring(s, 1, 1)

    if (s1 > 'a' && s1 < 'z')
        s <- paste(toupper(s1), substring(s, 2), sep='')

    s
}
