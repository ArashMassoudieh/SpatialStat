#include "environment.h"
#include "interface.h"
#include "grid.h"
#include "Distribution.h"
#include <QTableWidget>
#include <QCoreApplication>



Environment::Environment(QTableWidget *_outputwindow)
{
    outputwindow = _outputwindow;
}


bool Environment::Execute(const Command &cmd)
{
    if (outputwindow)
    {
        QTableWidgetItem *newItem = new QTableWidgetItem(QString::fromStdString(cmd.object_name + "." + cmd.command));
        outputwindow->insertRow(outputwindow->rowCount());
        outputwindow->setItem(outputwindow->rowCount()-1, 0, newItem);
        QProgressBar *pgbar = new QProgressBar();
        outputwindow->setCellWidget(outputwindow->rowCount()-1, 2,pgbar);
        current_progress_bar = pgbar;
        outputwindow->update();
        QCoreApplication::processEvents();
    }
    if (Grid::HasCommand(cmd.command))
    {
        if (cmd.Command_Structures[cmd.command].CommandType == command_type::creator)
        {
            Objects[cmd.object_name]=new Grid();
            Objects[cmd.object_name]->parent = this;
            dynamic_cast<Grid*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
            current_progress_bar->setValue(100);
            outputwindow->update();
            QCoreApplication::processEvents();
            return true;
        }
        if (cmd.Command_Structures[cmd.command].CommandType == command_type::modifier)
        {
            dynamic_cast<Grid*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
            current_progress_bar->setValue(100);
            outputwindow->update();
            QCoreApplication::processEvents();
            return true;
        }
    }
    if (CDistribution::HasCommand(cmd.command))
    {
        if (cmd.Command_Structures[cmd.command].CommandType == command_type::creator)
        {
            Objects[cmd.object_name]=new CDistribution();
            Objects[cmd.object_name]->parent = this;
            dynamic_cast<CDistribution*>(Objects[cmd.object_name])->Execute(cmd.command,cmd.arguments);
            current_progress_bar->setValue(100);
            outputwindow->update();
            QCoreApplication::processEvents();
            return true;
        }
        if (cmd.Command_Structures[cmd.command].CommandType == command_type::modifier)
        {
            dynamic_cast<CDistribution*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
            current_progress_bar->setValue(100);
            outputwindow->update();
            QCoreApplication::processEvents();
            return true;
        }
    }

    return false;
}
