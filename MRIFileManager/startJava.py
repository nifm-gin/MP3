import subprocess
option1 = '[ExportNifti] //tmp//Nifti'
option2 = '[ExportToMIA] PatientName-StudyName-CreationDate-SeqNumber-Protocol-SequenceName-AcquisitionTime'
option3 = 'CloseAfterExport'
option4 = 'NoLogExport'
option5 = '[LookAndFeel] com.sun.java.swing.plaf.windows.WindowsLookAndFeel'
option6 = 'NoExitSystem'
option7 = '[ExportOptions] 00013'
option8 = '[ProjectsDir] //home//omontigo//Documents//'

codexit=subprocess.call(['java','-Xms512m','-Xmx4096m','-jar','MRIManager.jar',option1,option2,option8,option3,option4,option7])
# subprocess.call('MRIManager.bat')
print(codexit)
