
Set objWshShell = WScript.CreateObject("WScript.Shell")
Set objFSO = WScript.CreateObject("Scripting.FileSystemObject")

strFolder = objFSO.getParentFolderName(WScript.ScriptFullName)

strText = strFolder & "\junction_list.txt"

Set objText = objFSO.OpenTextFile(strText, 2, true)

WScript.Echo "Junctions"

ListJunction objFSO.GetFolder("C:\cygwin_wm\")

objText.Close

Sub ListJunction(ByVal objMainFolder)
    Dim objSubFolder
    Dim objFile

    '// フォルダがあれば再帰
    For Each objSubFolder In objMainFolder.SubFolders
        ListJunction objSubFolder
    Next

    For Each objFile In objMainFolder.files
        strCmdLine = "cmd /c dir " & objFile.Path 
        Set objExecCmd = objWshShell.Exec(strCmdLine)

        strOut = objExecCmd.StdOut.ReadAll

        If InStr(strOut, "<JUNCTION>") Then
          WScript.Echo objFile.Path 
          strFile = Replace(objFile.Path, "\", "/")
          objText.WriteLine(strFile)
        End If
    Next
End Sub



