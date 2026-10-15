# Resolve an argument that has been renamed. `new` is the value of the current
# argument (NULL when it was not supplied) and `old` that of its deprecated alias.
# Supplying both is an error; the alias alone works with a warning.
.batchmixRenamedArg <- function(new, old, new_name, old_name, fn, required = FALSE) {
  if (!is.null(old)) {
    if (!is.null(new)) {
      stop(
        sprintf(
          "`%s` and its deprecated alias `%s` were both given to `%s()`; use `%s` only.",
          new_name, old_name, fn, new_name
        ),
        call. = FALSE
      )
    }
    warning(
      sprintf("The `%s` argument of `%s()` is deprecated; use `%s` instead.", old_name, fn, new_name),
      call. = FALSE
    )
    new <- old
  }
  if (required && is.null(new)) {
    stop(sprintf("argument \"%s\" is missing, with no default", new_name), call. = FALSE)
  }
  new
}
