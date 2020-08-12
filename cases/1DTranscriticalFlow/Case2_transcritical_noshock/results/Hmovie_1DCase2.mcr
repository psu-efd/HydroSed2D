#!MC 1000
$!VarSet |MFBD| = '/home/liu19/research/SWESed/cases/results'
$!VarSet |t| = 0
$!VarSet |dt| = 1
$!VarSet |realdt| = 1
$!VarSet |n| = 100
$!VarSet |time| = 0

# |StateAppend| = 0  --> if file exists, delete it first
#                 1 -->  if file exists, append
$!VarSet |StateAppend| = 0
#set the animation file format and name
$!EXPORTSETUP IMAGEWIDTH = 800
  EXPORTFORMAT = AVI
  ANIMATIONSPEED = 5
  QUALITY        = 100
  EXPORTFNAME = '|MFBD|/OneDCase2.avi'

$!LOOP |n|

  $!DRAWGRAPHICS FALSE
  $!OPENLAYOUT  "|MFBD|/oneDTestCase2.lay"
     ALTDATALOADINSTRUCTIONS = '"|MFBD|/h|t|tec.dat"'
     
  $!ATTACHTEXT 
  ANCHORPOS
    {
    X = 18
    Y = 86
    }
  COLOR = BLUE
  TEXT = 'Time = |Time| s' 
   
     
  $!DRAWGRAPHICS TRUE
  $!REDRAWALL

  $!IF |LOOP| == 1
    $!EXPORTSTART
    $!IF |StateAppend| == 0
      $!EXPORT
        APPEND = NO
    $!ENDIF
    $!IF |StateAppend| == 1
      $!EXPORT
        APPEND = YES
    $!ENDIF
  $!ENDIF
  $!IF |LOOP| > 1
    $!EXPORT
       APPEND = YES
  $!ENDIF
  $!VarSet |t| += |dt|
  $!VarSet |time| += |realdt|
$!ENDLOOP

$!EXPORTFINISH
$!RemoveVar |t|
$!RemoveVar |dt|
$!RemoveVar |realdt|
$!RemoveVar |StateAppend|
$!RemoveVar |n|
$!RemoveVar |MFBD|
