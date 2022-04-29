#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->pushButtonTest,SIGNAL(clicked()),this,SLOT(on_test()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_test()
{
    Command cmd("grid=CreateGrid(nx:100,ny:100,dx:0.01,dy:0.01)");
    environment.Execute(cmd);
}
