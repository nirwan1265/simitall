test_that("sequencing platform registry has the supported Illumina and PacBio families", {
  profiles <- simitall_sequencing_profiles()
  expect_true(all(c(
    "MiSeq", "MiSeq i100", "NextSeq 1000/2000", "NovaSeq X",
    "Sequel IIe HiFi", "Revio HiFi", "Vega HiFi"
  ) %in% profiles$platform))
})

test_that("platform names resolve to transparent simulator proxies", {
  miseq <- simitall:::.simitall_illumina_profile("MiSeq")
  nextseq <- simitall:::.simitall_illumina_profile("NextSeq 2000")
  revio <- simitall:::.simitall_pacbio_profile("Revio HiFi")
  expect_identical(miseq$art_system, "MSv3")
  expect_identical(nextseq$art_system, "NS50")
  expect_identical(revio$type, "HIFI")
  expect_error(simitall:::.simitall_illumina_profile("PromethION"), "Unsupported sequencing platform")
})

test_that("Illumina read-length bounds are checked before ART is required", {
  fasta <- tempfile(fileext = ".fa")
  writeLines(c(">chr1", "ACGTACGTACGT"), fasta)
  expect_error(
    sim_illumina_art(fasta, tempfile("miseq_"), cov = 10, platform = "NovaSeq X", readlen = 300),
    "150 bp"
  )
})
