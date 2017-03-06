{ file3d = $NF
  split(file3d, arr, ".")
  filexml = arr[1] ".xml"
  print "<scene> <model file='" file3d "'/> </scene>" > filexml
}
