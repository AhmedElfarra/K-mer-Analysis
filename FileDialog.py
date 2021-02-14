import sys
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog,QPushButton, QVBoxLayout, QDialog
from PyQt5.QtGui import QIcon


class App(QWidget):

    def __init__(self):
        super().__init__()
        # self.title = 'Kmer Counting'
        # self.left = 10
        # self.top = 10
        # self.width = 640
        # self.height = 480
        # self.edit = QLineEdit("")

        # self.button = QPushButton("Open file")
        # self.button = QPushButton("Open Directory")
        # # Create layout and add widgets
        # layout = QVBoxLayout()
        # layout.addWidget(self.edit)

        # layout.addWidget(self.button)
        # # Set dialog layout
        # self.setLayout(layout)
        # # Add button signal to greetings slot
        # self.button.clicked.connect(self.openFileNameDialog)
        # # self.edit = QLineEdit(self.openFileNameDialog())
        # self.initUI()
        
    # def initUI(self):
    #     self.setWindowTitle(self.title)
    #     self.setGeometry(self.left, self.top, self.width, self.height)
    
    def openFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            return fileName
        
    def openDirectory(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        dir_ = QFileDialog.getExistingDirectory(None, 'Select a folder:', 'C:\\', QFileDialog.ShowDirsOnly)
        
    def openFileNamesDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        files, _ = QFileDialog.getOpenFileNames(self,"QFileDialog.getOpenFileNames()", "","All Files (*);;Python Files (*.py)", options=options)
        if files:
            return files
    
    def saveFileDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;Text Files (*.txt)", options=options)
        if fileName:
            return fileName

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    ex.show()
    sys.exit(app.exec_())