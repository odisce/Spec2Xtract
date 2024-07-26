test_that(
  "load_files", {
    tempA <- rawrr::sampleFilePath()[1]
    tempB <- rep(tempA, 4)
    mapply(
      function(x) {
        temp <- load_files(x)
        expect_true(is.data.table(temp))
        expect_true(nrow(temp) == length(x))
        expect_true(
          all(
            c(
              "FileIndex",
              "FileExist",
              "file_name",
              "file_path"
            ) %in% names(temp)
          )
        )
        expect_true(all(seq_len(length(x)) %in% temp[, FileIndex]))
        expect_true("FileIndex" %in% names(temp))
        expect_true(all(isTRUE(temp[, FileExist])))
      },
      c(
        tempA,
        tempB
      )
    )
  }
)
