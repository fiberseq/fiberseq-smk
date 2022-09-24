#!/bin/tcsh -efx
# author : sjn
# date : Fri Sep 23 12:29:05 PDT 2022

# This script assumes the qc graphics are in the same dir as the <html-out>

if ( $#argv < 3 ) then
  printf "Expect $0 <sample-name> <html-out> <pdf>*\n"
  exit -1
endif

set sample=$1
set html=$2
set pdfs = "$argv[3-]"

cat <<__HTML__ >! $html
<html>
  <head>
    <style>
      table, th, td {
        border:1px solid black;
      }
      img {
        border: 1px solid #ddd;
        border-radius: 4px;
        padding: 5px;
        width: 150px;
      }
      img:hover {
        box-shadow: 0 0 2px 1px rgba(0, 140, 186, 0.5);
      }
    </style>
  </head>
  <body>
    <table style="width:100%">
    <caption><h2>$sample</h2></caption>
      <tr>
__HTML__

foreach pdf ($pdfs)
  convert $pdf $pdf:r.png
  printf '        <th>'$pdf:t:r'</th>' >>! $html
end

printf '\n      </tr>' >>! $html
printf '\n      <tr>\n' >>! $html
foreach pdf ($pdfs)
  printf '        <td><a href="'$pdf:t'"><img src="'$pdf:t:r.png'"></img></a></td>' >>! $html
end
printf '\n' >>! $html

cat << __HTML__ >>! $html
      </tr>
    </table>
  </body>
</html>
__HTML__

exit 0
