Output docs/demo.gif

Set FontSize 16
Set Width 1200
Set Height 760
Set Padding 12
Set Framerate 24
Set Theme "Catppuccin Mocha"

Hide
Type "pixelfold 1UBQ --no-surface" Enter
Sleep 3s
Show

Sleep 1200ms
Type "wwww"   Sleep 700ms      # rotate (pitch)
Type "dddddd" Sleep 700ms      # rotate (yaw)
Type "++"     Sleep 700ms      # zoom in
Type "b"      Sleep 1400ms     # b-factor coloring
Type "b"      Sleep 600ms      # back to secondary-structure coloring
Type "2"      Sleep 1200ms     # C-alpha backbone
Type "1"      Sleep 600ms      # all atoms
Type "h"      Sleep 1400ms     # hydrogen bonds
Type "n"      Sleep 1800ms     # H-bond network overlay
Type "hn"     Sleep 400ms      # clear both
Type "v"      Sleep 2500ms     # solvent-accessible surface (computed on demand)
Type "q"
