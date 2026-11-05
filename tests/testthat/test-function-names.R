test_that("package functions use lowerCamelCase", {
    ns <- asNamespace("hpgR")
    functions <- Filter(function(name) is.function(get(name, ns, inherits = FALSE)),
                        ls(ns, all.names = TRUE))
    expect_true(all(grepl("^[a-z][A-Za-z0-9]*$", functions)))
    checkNames <- function(x) {
        if (!is.call(x) && !is.expression(x) && !is.pairlist(x))
            return(invisible(NULL))
        if (is.call(x) && is.symbol(x[[1L]]) &&
            as.character(x[[1L]]) %in% c("<-", "=") && length(x) == 3L &&
            is.symbol(x[[2L]]) && is.call(x[[3L]]) &&
            identical(x[[3L]][[1L]], as.name("function")))
            expect_match(as.character(x[[2L]]), "^[a-z][A-Za-z0-9]*$")
        for (i in seq_along(x)) {
            if (identical(x[[i]], quote(expr = ))) next
            if (is.call(x[[i]]) || is.expression(x[[i]]) || is.pairlist(x[[i]]))
                checkNames(x[[i]])
        }
        invisible(NULL)
    }
    for (name in functions) checkNames(body(get(name, ns, inherits = FALSE)))
})
