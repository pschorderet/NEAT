The officially sanctioned way is now the iconutil command instead.

From Apple's guidelines:

After you’ve created the necessary app icon assets, place them in folder a named icon.iconset. To create an .icns file, use iconutil in Terminal. Terminal is located in /Applications/Utilities/Terminal. Enter the command iconutil -c icns <iconset filename>, where <iconset filename> is the path to the .iconset folder. You must use iconutil, not Icon Composer, to create high-resolution .icns files.
There's another relevant Apple doc that goes more in depth: High Resolution Resources.

For reference, the complete set of icons:

icon_16x16.png
icon_16x16@2x.png
icon_32x32.png
icon_32x32@2x.png
icon_128x128.png
icon_128x128@2x.png
icon_256x256.png
icon_256x256@2x.png
icon_512x512.png
icon_512x512@2x.png