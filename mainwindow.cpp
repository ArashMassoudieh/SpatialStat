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
    Command cmd = Command("grid=CreateGrid(nx:100,ny:100,dx:0.01,dy:0.01)");
    //environment.Execute(cmd);
    cmd = Command("dist=CreateDistribution(type:lognormal,p0:0,p1:1)");
    //environment.Execute(cmd1);
    cmd = Command("dist*WriteToFile(filename:test.txt,nbins:50");
    //environment.Execute(cmd2);
}
