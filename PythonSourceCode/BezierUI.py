from PySide2.QtGui import *
from PySide2.QtCore import  *
from PySide2.QtWidgets import *
import cmds_weight as weight

def q_add(layout,*elements):
    for element in elements:
        if isinstance(element,QLayout):
            layout.addLayout(element)
        elif isinstance(element,QWidget):
            layout.addWidget(element)
    return layout

def q_button(text,action):
    button = QPushButton(text)
    button.clicked.connect(action)
    return button


class BezierWindow(QWidget):
    valueChanged = Signal()

    def __init__(self):
        QWidget.__init__(self)
        self.points = [[0.0,1.0],[0.5,1.0],[0.5,0.0],[1.0,0.0]]
        self.__movePoint = 0
        self.setFixedSize(500,500)
        self.__mirror = False
        self.__magnet = False
        self.setWindowTitle("Bezier Curve")
        self.valueChanged.connect(self.getParameterOfBezier)

    def paintEvent(self,event):
        QWidget.paintEvent(self,event)
        # paint background
        painter = QPainter(self)
        painter.setBrush(QBrush(QColor(120,120,120),Qt.SolidPattern))
        painter.setPen(QPen(QColor(0,0,0),1,Qt.SolidLine))
        painter.drawRect(0,0,self.width()-1,self.height()-1)

        # paint curve
        painter.setBrush(QBrush(QColor(100,100,100),Qt.SolidPattern))
        points = [QPointF((self.width()-1)*p[0],(self.height()-1)*p[1]) for p in self.points]
        path = QPainterPath()
        path.moveTo(self.width()-1,self.height()-1)
        path.lineTo(points[0])
        path.cubicTo(*points[1:])
        painter.drawPath(path)

        # paint grid
        painter.setPen(QPen(QColor(200,200,200),1,Qt.DotLine))
        width_step = (self.width()-1)/10
        height_step = (self.height()-1)/10
        for i in range(1,10):
            width = width_step*i
            height = height_step*i
            painter.drawLine(width,0,width,self.height())
            painter.drawLine(0,height,self.width(),height)

        # paint control points
        painter.setPen(QPen(QColor(0, 0, 0), 1, Qt.SolidLine))
        painter.setBrush(QBrush(QColor(200, 200, 200), Qt.SolidPattern))
        painter.drawEllipse(points[1], 8, 8)
        painter.drawEllipse(points[2], 8, 8)
        # paint control lines
        painter.setPen(QPen(QColor(200,200,200),1,Qt.DashLine))
        painter.drawLine(points[0],points[1])
        painter.drawLine(points[3],points[2])
        # paint axis values
        painter.setPen(QPen(QColor(0,0,0),2,Qt.SolidLine))
        for i in range(1, 11):
            x_value = i * 0.1 * (self.width() - 1)
            painter.drawText(x_value-15, self.height()-2, f"{i * 0.1:.1f}")
        for i in range(1, 11):
            y_value = (1-i * 0.1) * (self.height() - 1)
            painter.drawText(3, y_value+10, f"{i * 0.1:.1f}")

        painter.end()

    # use for getting the selected control points
    def mousePressEvent(self,event):
        self.setFocus()
        QWidget.mousePressEvent(self,event)
        points = [QPointF((self.width()-1)*p[0],(self.height()-1)*p[1]) for p in self.points]
        for i in [1,2]:
            p = QPointF(event.pos()) - points[i]
            length = (p.x() ** 2 + p.y() ** 2) ** 0.5
            if length < 10:
                self.__movePoint =i
                return
        self.__movePoint =0

    def keyPressEvent(self, event):
        QWidget.keyPressEvent(self, event)
        if event.key() == Qt.Key_X:
            self.__magnet = True
        if event.modifiers() == Qt.ControlModifier:
            self.__mirror = True

    def keyReleaseEvent(self, event):
        QWidget.keyReleaseEvent(self, event)
        self.__mirror = False
        self.__magnet = False

    def mouseMoveEvent(self,event):
        QWidget.mouseMoveEvent(self,event)
        if self.__movePoint ==1:
            p = QPointF(event.pos())
            x = max(min(float(p.x()/(self.width()-1)),1.0 ),0.0)
            y = max(min(float(p.y()/(self.height()-1)),1.0 ),0.0)
            if self.__magnet:
                x = round(x*10)/10.0
                y = round(y*10)/10.0
            if self.__mirror:
                mx = (1-x)
                my = (1-y)
                self.points[2] = [mx,my]
            self.points[1]=[x,y]
            self.update()
            self.valueChanged.emit()

        elif self.__movePoint ==2:
            p = QPointF(event.pos())
            x = max(min(float(p.x()/(self.width()-1)),1.0 ),0.0)
            y = max(min(float(p.y()/(self.height()-1)),1.0 ),0.0)
            if self.__magnet:
                x = round(x*10)/10.0
                y = round(y*10)/10.0
            if self.__mirror:
                mx = (1-x)
                my = (1-y)
                self.points[1] = [mx,my]
            self.points[2]=[x,y]
            self.update()
            self.valueChanged.emit()

    def getParameterOfBezier(self):
        Px = [point[0] for point in self.points]
        Py = [1-point[1] for point in self.points]
        return [Px,Py]

class FloatSliderGroup(QHBoxLayout):
    valueChanged = Signal(float)

    def __init__(self):
        QHBoxLayout.__init__(self)
        self.slider = QSlider(Qt.Horizontal)
        self.spin = QDoubleSpinBox()
        self.addWidget(self.spin)
        self.addWidget(self.slider)
        self.set_range(0, 1)
        self.spin.valueChanged.connect(self.convert_slider)
        self.slider.valueChanged.connect(self.convert_spin)
        self.spin.setSingleStep(0.01)
        self.spin.setDecimals(3)

    def convert_spin(self, value):
        self.spin.setValue(value/1000.0)
        self.valueChanged.emit(self.spin.value())

    def convert_slider(self, value):
        self.slider.setValue(int(round(1000*value)))

    def set_range(self, min_value, max_value):
        self.spin.setRange(min_value, max_value)
        self.slider.setRange(min_value*1000, max_value*1000)

    def value(self):
        return self.spin.value()

    def set_value(self, value):
        self.spin.setValue(value)


class BezierToolUI(QDialog):
    # main UI
    def __init__(self,parent):
        # base set up
        QDialog.__init__(self, parent)
        self.setWindowTitle("Dante's Bezier Weight Tool")
        self.resize(500, 800)
        maninLayout=QVBoxLayout()
        self.setLayout(maninLayout)
        # bezier window
        self.bezierWindow = BezierWindow()
        maninLayout.addWidget(self.bezierWindow)
        # add stretch for the rest part
        maninLayout.addStretch()
        H_layout = QHBoxLayout()
        maninLayout.addLayout(H_layout)
        H_layout.addStretch()
        # prepare a specific layout
        Args_form_layout = QFormLayout()
        H_layout.addLayout(Args_form_layout)
        H_layout.addStretch()

        # first real time button
        self.mode = QCheckBox()
        Args_form_layout.addRow(U"Real Time", self.mode)

        # option buttons
        self.type = QButtonGroup()
        Grid_layout = QGridLayout()
        for i, text in enumerate([u"Single Joint", u"Soft", u"Joints Chain"]):
            element = QRadioButton(text)
            self.type.addButton(element,i)
            Grid_layout.addWidget(element,i/3,i%3,1,1)
        self.type.button(0).setChecked(True)
        Args_form_layout.addRow(u"Type",Grid_layout)

        # axis buttons
        self.Axis = QButtonGroup()
        Grid_layout = QGridLayout()
        for i, text in enumerate([u"X", u"Y", u"Z",u"distance"]):
            element = QRadioButton(text)
            self.Axis.addButton(element, i)
            Grid_layout.addWidget(element, i / 3, i % 3, 1, 1)
        self.Axis.button(0).setChecked(True)
        Args_form_layout.addRow(u"Axis", Grid_layout)

        # radius slider
        self.radius = FloatSliderGroup()
        self.radius.set_range(0,100)
        self.radius.set_value(3)
        Args_form_layout.addRow(u"Radius",self.radius)
        maninLayout.addStretch()

        # calculate button
        calculate_button =q_button("calculate the weight",self.calculate)
        maninLayout.addWidget(calculate_button)

        self.singleInfoList = []
        self.chainInfoList = []

        self.type.buttonClicked.connect(self.get_Args)
        self.Axis.buttonClicked.connect(self.get_Args)
        self.bezierWindow.valueChanged.connect(self.realTime_Calculate)
        self.radius.valueChanged.connect(self.realTime_Calculate)


    # get the args for compute
    def get_Args(self):
        bezier_list = self.bezierWindow.getParameterOfBezier()
        Px = bezier_list[0]
        Py = bezier_list[1]
        radius = self.radius.value()
        type = self.type.id(self.type.checkedButton())
        axis = self.Axis.id(self.Axis.checkedButton())
        return dict(
            Px=Px,
            Py=Py,
            radius = radius,
            type =type,
            axis = axis
        )

    def calculate(self):
        if self.mode.isChecked():
            self.realTime_Calculate()
        else:
            Args = self.get_Args()
            type = Args["type"]
            Px = Args["Px"]
            Py = Args["Py"]
            radius = Args["radius"]
            axis = Args["axis"]
            if type==0:
                self.singleInfoList = weight.get_single_joint_info()
                weight.paint_SingleWeights(px=Px,py=Py,radius=radius,axis=axis,infoList=self.singleInfoList)
            if type==1:
                weight.paint_soft(px=Px,py=Py,radius=None)
            if type==2:
                self.chainInfoList = weight.get_chain_info()
                weight.paint_Chain_weights(px=Px,py=Py,radius=radius,axis=axis,ChainInfoList=self.chainInfoList)


    def realTime_Calculate(self):
        if self.mode.isChecked():
            Args = self.get_Args()
            type = Args["type"]
            Px = Args["Px"]
            Py = Args["Py"]
            radius = Args["radius"]
            axis = Args["axis"]
            if type == 0:
                weight.paint_SingleWeights(px=Px, py=Py, radius=radius, axis=axis,infoList=self.singleInfoList)
            if type == 1:
                weight.setToVertexMode()
                weight.paint_soft(px=Px, py=Py, radius=radius)
            if type == 2:
                weight.paint_Chain_weights(px=Px, py=Py, radius=radius, axis=axis,ChainInfoList=self.chainInfoList)




def get_Top():
    top = QApplication.activeWindow()
    if top is None:
        return
    while True:
        parent = top.parent()
        if parent is None:
            return top
        top = parent

    # Function to show the dialog window
def show():
    global window
    if window is None:
        # Create a new QDialog window with the top-level application window as its parent
        window = BezierToolUI(parent=get_Top())
    # Show the dialog window
    window.show()

if __name__ == "__main__":
    my_app = QApplication([])
    show()
    my_app.exec_()

window = None