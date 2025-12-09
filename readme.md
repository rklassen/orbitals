
# pip-installing python libraries in blender

Option 1: Disable SIP temporarily (requires restart)
This is the most reliable but requires rebooting:

Restart your Mac and hold Cmd + R during boot (Intel) or hold power button and select Options (Apple Silicon)
In Recovery Mode, open Terminal from Utilities menu
Run: csrutil disable
Restart normally
Run the copy commands
Re-enable SIP by repeating steps 1-2 and running: csrutil enable


cd /Applications/Blender.app/Contents/Resources/5.0/python/bin

# Install to temp location first
./python3.11 -m pip install scipy --target ~/scipy_temp

# Then copy it over
sudo cp -R ~/scipy_temp/scipy /Applications/Blender.app/Contents/Resources/5.0/python/lib/python3.11/site-packages/
sudo cp -R ~/scipy_temp/scipy*.dist-info /Applications/Blender.app/Contents/Resources/5.0/python/lib/python3.11/site-packages/

# Fix ownership
sudo chown -R $(whoami):staff /Applications/Blender.app/Contents/Resources/5.0/python/lib/python3.11/site-packages/scipy*

# Verify it's there
ls -l /Applications/Blender.app/Contents/Resources/5.0/python/lib/python3.11/site-packages/ | grep scipy

# Clean up temp folder
rm -rf ~/scipy_temp