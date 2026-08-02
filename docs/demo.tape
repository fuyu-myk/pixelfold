Output docs/demo.gif

Set FontSize 11
Set Width 1280
Set Height 800
Set Padding 12
Set Framerate 24
Set Theme "Catppuccin Mocha"

Hide
Type "export PATH=$PWD/target/release:$PATH" Enter
Type "clear" Enter
Show

Type "pixelfold 4DFR --protocol half-block" Enter
Sleep 4s

Sleep 800ms
Type "wwww"   Sleep 700ms      # rotate (pitch)
Type "dddddd" Sleep 900ms      # rotate (yaw)
Type "++"     Sleep 700ms      # zoom in
Type "b"      Sleep 1200ms     # b-factor colouring
Type "b"      Sleep 500ms      # back to secondary structure
Type "2"      Sleep 1000ms     # C-alpha backbone
Type "1"      Sleep 500ms      # all atoms
Type "h"      Sleep 1200ms     # hydrogen bonds
Type "h"      Sleep 400ms      # clear them

Type "g"      Sleep 2500ms     # linked network pane
Type "e"      Sleep 1600ms     # colour nodes by burial (RSA)
Type "l"      Sleep 1800ms     # spatial layout, tracks the camera
Type "dddd"   Sleep 1200ms     # rotate: the graph rotates with the structure
Type "l"      Sleep 1000ms     # back to the force-directed layout
Type "y"      Sleep 1400ms     # chain filter: intra-chain only
Type "y"      Sleep 1400ms     # inter-chain interface only
Type "y"      Sleep 800ms      # all edges again

Type "q"
