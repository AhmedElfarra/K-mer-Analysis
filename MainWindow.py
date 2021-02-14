import sys
import os
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5 import QtCore
from PyQt5.QtWidgets import (
    QMainWindow,
    QProgressBar,
    QApplication,
    QHBoxLayout,
    QPushButton,
    QWidget,
    QLineEdit,
    QDialog,
    QGridLayout,
    QFileDialog,
    QLabel
)

class Window(QWidget):
    def __init__(self):
        super().__init__()
        self.title = "K-mer Analysis Tool"
        self.left = 500 
        self.top = 200 
        self.width = 550
        self.height = 420
        self.iconName = "pngegg.png"
        self.setWindowIcon(QIcon("pngegg.png"))
        
        self.InitWindow()
        
    def InitWindow(self):
        self.setWindowTitle(self.title)
        self.setWindowIcon(QIcon("self.iconName"))
        self.setGeometry(self.left,self.top,self.width,self.height)
        self.setFixedSize(self.width,self.height)
        self.CreateLayout()
        
        self.show()
        
    def CreateLayout(self):
        self.kLabel = QLabel('Enter The Value For K:')
        self.folder1 = QPushButton("Open Folder 1:")
        self.folder2 = QPushButton("Open Folder 2:")
        self.kmerButton = QPushButton("OK")
        self.processData = QPushButton("Process Data")
        self.download = QPushButton("Download")
     
        
        kLabelEdit = QLineEdit()
        folderEdit1 = QLineEdit()
        folderEdit2 = QLineEdit()
        
        layout = QGridLayout()
        layout.setSpacing(10)
        labelImage = QLabel(self)
        pixmap = QPixmap("pngegg.png")
        pixmap5 = pixmap.scaled(240, 240)
        labelImage.setPixmap(pixmap5)
        layout.addWidget(labelImage,0,1)
        
        layout.addWidget(self.kLabel, 1, 0)
        layout.addWidget(kLabelEdit, 1, 1)
        layout.addWidget(self.kmerButton, 1, 2)
        
        
        layout.addWidget(self.folder1, 2, 0)
        self.folder1.clicked.connect(self.FileDialog)
        layout.addWidget(folderEdit1, 2, 1) 

        layout.addWidget(self.folder2, 3, 0)
        layout.addWidget(folderEdit2, 3, 1)
        self.folder2.clicked.connect(self.FileDialog)
                
        layout.addWidget(self.processData, 6, 1)
        layout.addWidget(self.download, 6, 2)

        self.setLayout(layout)

            
    def FileDialog(self,directory='', forOpen=True, fmt='', isFolder=True):
        ROOT_DIR = os.path.dirname(os.path.abspath("top_level_file.txt"))
        
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        options |= QFileDialog.DontUseCustomDirectoryIcons
        dialog = QFileDialog()
        dialog.setOptions(options)

        dialog.setFilter(dialog.filter() | QtCore.QDir.Hidden)

        # ARE WE TALKING ABOUT FILES OR FOLDERS
        if isFolder:
            dialog.setFileMode(QFileDialog.DirectoryOnly)
        else:
            dialog.setFileMode(QFileDialog.AnyFile)
        # OPENING OR SAVING
        dialog.setAcceptMode(QFileDialog.AcceptOpen) if forOpen else dialog.setAcceptMode(QFileDialog.AcceptSave)

        # SET FORMAT, IF SPECIFIED
        if fmt != '' and isFolder is False:
            dialog.setDefaultSuffix(fmt)
            dialog.setNameFilters([f'{fmt} (*.{fmt})'])

        # SET THE STARTING DIRECTORY
        if directory != '':
            dialog.setDirectory(str(directory))
        else:
            dialog.setDirectory(str(ROOT_DIR))
        

        if dialog.exec_() == QDialog.Accepted:
            path = dialog.selectedFiles()[0]  # returns a list
            return path
        else:
            return ''
        

    def valueK(self):
        return self.kCount.text()
        
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = Window()
    sys.exit(app.exec_())