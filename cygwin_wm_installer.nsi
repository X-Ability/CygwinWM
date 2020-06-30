; NSIS Script

!include MUI2.nsh

!verbose 4

; /SOLID������ƃG���[�ɂȂ�(�T�C�Y���傫�����邩��H)
SetCompressor lzma

; �C���X�g�[���[�̖��O
Name "Cygwin for Winmostar"

; �C���X�g�[���[�̃t�@�C����
OutFile "cygwin_wm.exe"

; �C���X�g�[����f�B���N�g��
InstallDir "C:\cygwin_wm"

; �K�v�Ȏ��s�����iuser/admin�j
RequestExecutionLevel user

; �C���X�g�[�����ɕ\��������̂��w��
; �㏑���֎~
!define MUI_PAGE_CUSTOMFUNCTION_LEAVE ValidateDirectory
!insertmacro MUI_PAGE_DIRECTORY
!define MUI_PAGE_CUSTOMFUNCTION_LEAVE CreateSymlinks
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_LANGUAGE "English"

; �C���X�g�[�����̋���
Section "Install" ; ���̖��O�͏d�v�ł͂Ȃ�

  ; �L�����Z���ł���悤��
  GetDlgItem $0 $HWNDPARENT 2
  EnableWindow $0 1

  SectionIn RO

  ; �C���X�g�[����f�B���N�g�����w��
  SetOutPath $INSTDIR

  SetOverwrite off
  File /a /r "C:\cygwin_wm\*.*"

SectionEnd

Function ValidateDirectory
  StrCpy $0 $INSTDIR 16
  StrCmp $0 "C:\Program Files" folderwarning
  StrCpy $0 $INSTDIR 22
  StrCmp $0 "C:\Program Files (x86)" folderwarning
  Goto folderok

  folderwarning:

  MessageBox MB_OK|MB_ICONSTOP "Cygwin does not work properly with folders that require administrator privileges for writing."
  Abort

  folderok:

  ${If} ${FileExists} `$INSTDIR\*.*`
    MessageBox MB_OK|MB_ICONSTOP "Installation folder already exists. Move or delete the folder before installation"
    Abort
  ${EndIf}

FunctionEnd

Function CreateSymlinks
  ExecWait '"$INSTDIR\make_symlink.cmd"'
FunctionEnd

