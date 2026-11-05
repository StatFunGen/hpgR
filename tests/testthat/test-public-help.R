test_that('each exported function has a signature, arguments and return contract', {
    pages <- tools::Rd_db('hpgR', lib.loc = .libPaths())
    tag <- function(page, name) Filter(function(x)
        identical(attr(x,'Rd_tag'),paste0('\\',name)),page)
    content <- function(x) paste(unlist(x),collapse='')
    for (name in getNamespaceExports('hpgR')) {
        matching <- Filter(function(p) name %in%
            vapply(tag(p,'alias'),content,''),pages)
        expect_equal(length(matching),1L,info=name)
        page <- matching[[1L]]
        expect_true(length(tag(page,'usage'))>0L,info=name)
        expect_true(length(tag(page,'value'))>0L,info=name)
        calls <- parse(text=content(tag(page,'usage')))
        usage_names <- vapply(calls,function(x) as.character(x[[1L]]),'')
        expect_true(name %in% usage_names,info=name)
        params <- names(formals(getExportedValue('hpgR',name)))
        if (length(params)) {
            entries <- tag(page,'arguments')[[1L]]
            items <- Filter(function(x) identical(attr(x,'Rd_tag'),'\\item'),entries)
            documented <- unlist(strsplit(vapply(items,function(x)content(x[[1L]]),''),
                                          ',[[:space:]]*'))
            expect_true(all(params %in% documented),info=name)
        }
    }
})
