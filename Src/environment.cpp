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
    FunctionOutPut output;
    if (Grid::HasCommand(cmd.command))
    {
        if (Objects.count(cmd.object_name)==0)
        {
            Objects[cmd.object_name] = new Grid();
            Objects[cmd.object_name]->parent = this;
        }
        if (cmd.Command_Structures[cmd.command].Output==object_type::grid)
        {
            output = dynamic_cast<Grid*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
            Objects[cmd.object_name] = output.output;
        }
        else if (cmd.Command_Structures[cmd.command].Output==object_type::distribution)
        {
            output = dynamic_cast<Grid*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
            Objects[cmd.object_name] = output.output;
        }
        else if (cmd.Command_Structures[cmd.command].Output==object_type::timeseries)
        {
            output = dynamic_cast<Grid*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
            Objects[cmd.object_name] = output.output;
        }
        else
        {
            output = dynamic_cast<Grid*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
        }
        Objects[cmd.object_name]->parent = this;
        current_progress_bar->setValue(100);
        outputwindow->update();
        QCoreApplication::processEvents();
        return output.success;
    }
    if (CDistribution::HasCommand(cmd.command))
    {
        if (Objects.count(cmd.object_name)==0)
        {
            Objects[cmd.object_name] = new CDistribution();
            Objects[cmd.object_name]->parent = this;
        }
        if (cmd.Command_Structures[cmd.command].Output==object_type::grid)
        {
            output = dynamic_cast<CDistribution*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
            Objects[cmd.object_name] = output.output;
        }
        else if (cmd.Command_Structures[cmd.command].Output==object_type::distribution)
        {
            output = dynamic_cast<CDistribution*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
            Objects[cmd.object_name] = output.output;
        }
        else if (cmd.Command_Structures[cmd.command].Output==object_type::timeseries)
        {
            output = dynamic_cast<CDistribution*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
            Objects[cmd.object_name] = output.output;
         }
        else
        {
            output = dynamic_cast<CDistribution*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
        }
        Objects[cmd.object_name]->parent = this;
        current_progress_bar->setValue(100);
        outputwindow->update();
        QCoreApplication::processEvents();
        return output.success;
    }
    if (TimeSeriesD::HasCommand(cmd.command))
    {
        if (Objects.count(cmd.object_name)==0)
        {
            Objects[cmd.object_name] = new TimeSeriesD();
            Objects[cmd.object_name]->parent = this;
        }
        if (cmd.Command_Structures[cmd.command].Output==object_type::grid)
        {
            output = dynamic_cast<CDistribution*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
            Objects[cmd.object_name] = output.output;
        }
        else if (cmd.Command_Structures[cmd.command].Output==object_type::distribution)
        {
            output = dynamic_cast<CDistribution*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
            Objects[cmd.object_name] = output.output;
        }
        else if (cmd.Command_Structures[cmd.command].Output==object_type::timeseries)
        {
            output = dynamic_cast<CDistribution*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
            Objects[cmd.object_name] = output.output;
        }
        else
        {
            output = dynamic_cast<CDistribution*>(Objects[cmd.object_name])->Execute(cmd.command, cmd.arguments);
        }
        Objects[cmd.object_name]->parent = this;
        current_progress_bar->setValue(100);
        outputwindow->update();
        QCoreApplication::processEvents();
        return output.success;

    }

    return false;
}
