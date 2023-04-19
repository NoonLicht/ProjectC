#pragma once
#include <string.h>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <ctime>
#include "MySourceLib.cpp"
#include <msclr/marshal_cppstd.h>

#using <system.windows.forms.dll>
#using <Microsoft.VisualBasic.dll>

using namespace std;
using namespace msclr::interop;

namespace MyForm1 {

	using namespace System;
	using namespace System::IO;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Сводка для MyForm
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		MyForm(void)
		{
			InitializeComponent();
			//
			//TODO: добавьте код конструктора
			//
		}

	protected:
		/// <summary>
		/// Освободить все используемые ресурсы.
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::ToolStrip^ toolStrip1;
	private: System::Windows::Forms::ToolStripButton^ создатьToolStripButton;
	private: System::Windows::Forms::ToolStripButton^ открытьToolStripButton;
	private: System::Windows::Forms::ToolStripButton^ сохранитьToolStripButton;
	private: System::Windows::Forms::ToolStripButton^ печатьToolStripButton;
	private: System::Windows::Forms::ToolStripSeparator^ toolStripSeparator;
	private: System::Windows::Forms::ToolStripButton^ вырезатьToolStripButton;
	private: System::Windows::Forms::ToolStripButton^ копироватьToolStripButton;
	private: System::Windows::Forms::ToolStripButton^ вставкаToolStripButton;
	private: System::Windows::Forms::ToolStripSeparator^ toolStripSeparator1;
	private: System::Windows::Forms::ToolStripButton^ справкаToolStripButton;
	private: System::Windows::Forms::ToolStripButton^ toolStripButton1;
	private: System::Windows::Forms::ToolStripButton^ toolStripButton2;
	private: System::Windows::Forms::ToolStripButton^ toolStripButton3;
	private: System::Windows::Forms::ToolStripButton^ toolStripButton4;
	private: System::Windows::Forms::ToolStripButton^ toolStripButton5;
	private: System::Windows::Forms::ToolStripButton^ toolStripButton6;
	private: System::Windows::Forms::ToolStripButton^ toolStripButton7;
	private: System::Windows::Forms::ToolStripButton^ toolStripButton8;
	private: System::Windows::Forms::ToolStripButton^ toolStripButton9;
	private: System::Windows::Forms::RichTextBox^ richTextBox1;
	private: System::Windows::Forms::TreeView^ treeView1;
	private: System::Windows::Forms::DataVisualization::Charting::Chart^ chart1;
	protected:

	private:
		/// <summary>
		/// Обязательная переменная конструктора.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Требуемый метод для поддержки конструктора — не изменяйте 
		/// содержимое этого метода с помощью редактора кода.
		/// </summary>
		void InitializeComponent(void)
		{
			System::ComponentModel::ComponentResourceManager^ resources = (gcnew System::ComponentModel::ComponentResourceManager(MyForm::typeid));
			System::Windows::Forms::TreeNode^ treeNode1 = (gcnew System::Windows::Forms::TreeNode(L"Нормальный закон распределения"));
			System::Windows::Forms::TreeNode^ treeNode2 = (gcnew System::Windows::Forms::TreeNode(L"Распределение Вейбулла"));
			System::Windows::Forms::TreeNode^ treeNode3 = (gcnew System::Windows::Forms::TreeNode(L"Метод максимального правдоподобия", gcnew cli::array< System::Windows::Forms::TreeNode^  >(2) {
				treeNode1,
					treeNode2
			}));
			System::Windows::Forms::TreeNode^ treeNode4 = (gcnew System::Windows::Forms::TreeNode(L"Нормальный закон"));
			System::Windows::Forms::TreeNode^ treeNode5 = (gcnew System::Windows::Forms::TreeNode(L"Вейбулл"));
			System::Windows::Forms::TreeNode^ treeNode6 = (gcnew System::Windows::Forms::TreeNode(L"Метод наименьших квадратов", gcnew cli::array< System::Windows::Forms::TreeNode^  >(2) {
				treeNode4,
					treeNode5
			}));
			System::Windows::Forms::TreeNode^ treeNode7 = (gcnew System::Windows::Forms::TreeNode(L"Регрессивный анализ"));
			System::Windows::Forms::TreeNode^ treeNode8 = (gcnew System::Windows::Forms::TreeNode(L"Статистическое оценивание", gcnew cli::array< System::Windows::Forms::TreeNode^  >(3) {
				treeNode3,
					treeNode6, treeNode7
			}));
			System::Windows::Forms::TreeNode^ treeNode9 = (gcnew System::Windows::Forms::TreeNode(L"Обработка данных на ЭВМ", gcnew cli::array< System::Windows::Forms::TreeNode^  >(1) { treeNode8 }));
			System::Windows::Forms::DataVisualization::Charting::ChartArea^ chartArea1 = (gcnew System::Windows::Forms::DataVisualization::Charting::ChartArea());
			System::Windows::Forms::DataVisualization::Charting::Legend^ legend1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Legend());
			System::Windows::Forms::DataVisualization::Charting::Series^ series1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			this->toolStrip1 = (gcnew System::Windows::Forms::ToolStrip());
			this->создатьToolStripButton = (gcnew System::Windows::Forms::ToolStripButton());
			this->открытьToolStripButton = (gcnew System::Windows::Forms::ToolStripButton());
			this->сохранитьToolStripButton = (gcnew System::Windows::Forms::ToolStripButton());
			this->печатьToolStripButton = (gcnew System::Windows::Forms::ToolStripButton());
			this->toolStripSeparator = (gcnew System::Windows::Forms::ToolStripSeparator());
			this->вырезатьToolStripButton = (gcnew System::Windows::Forms::ToolStripButton());
			this->копироватьToolStripButton = (gcnew System::Windows::Forms::ToolStripButton());
			this->вставкаToolStripButton = (gcnew System::Windows::Forms::ToolStripButton());
			this->toolStripSeparator1 = (gcnew System::Windows::Forms::ToolStripSeparator());
			this->справкаToolStripButton = (gcnew System::Windows::Forms::ToolStripButton());
			this->toolStripButton1 = (gcnew System::Windows::Forms::ToolStripButton());
			this->toolStripButton2 = (gcnew System::Windows::Forms::ToolStripButton());
			this->toolStripButton3 = (gcnew System::Windows::Forms::ToolStripButton());
			this->toolStripButton4 = (gcnew System::Windows::Forms::ToolStripButton());
			this->toolStripButton5 = (gcnew System::Windows::Forms::ToolStripButton());
			this->toolStripButton6 = (gcnew System::Windows::Forms::ToolStripButton());
			this->toolStripButton7 = (gcnew System::Windows::Forms::ToolStripButton());
			this->toolStripButton8 = (gcnew System::Windows::Forms::ToolStripButton());
			this->toolStripButton9 = (gcnew System::Windows::Forms::ToolStripButton());
			this->richTextBox1 = (gcnew System::Windows::Forms::RichTextBox());
			this->treeView1 = (gcnew System::Windows::Forms::TreeView());
			this->chart1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Chart());
			this->toolStrip1->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->BeginInit();
			this->SuspendLayout();
			// 
			// toolStrip1
			// 
			this->toolStrip1->ImageScalingSize = System::Drawing::Size(20, 20);
			this->toolStrip1->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(19) {
				this->создатьToolStripButton,
					this->открытьToolStripButton, this->сохранитьToolStripButton, this->печатьToolStripButton, this->toolStripSeparator, this->вырезатьToolStripButton,
					this->копироватьToolStripButton, this->вставкаToolStripButton, this->toolStripSeparator1, this->справкаToolStripButton, this->toolStripButton1,
					this->toolStripButton2, this->toolStripButton3, this->toolStripButton4, this->toolStripButton5, this->toolStripButton6, this->toolStripButton7,
					this->toolStripButton8, this->toolStripButton9
			});
			this->toolStrip1->Location = System::Drawing::Point(0, 0);
			this->toolStrip1->Name = L"toolStrip1";
			this->toolStrip1->Size = System::Drawing::Size(1478, 27);
			this->toolStrip1->TabIndex = 0;
			this->toolStrip1->Text = L"toolStrip1";
			// 
			// создатьToolStripButton
			// 
			this->создатьToolStripButton->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->создатьToolStripButton->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"создатьToolStripButton.Image")));
			this->создатьToolStripButton->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->создатьToolStripButton->Name = L"создатьToolStripButton";
			this->создатьToolStripButton->Size = System::Drawing::Size(29, 24);
			this->создатьToolStripButton->Text = L"&Создать";
			// 
			// открытьToolStripButton
			// 
			this->открытьToolStripButton->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->открытьToolStripButton->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"открытьToolStripButton.Image")));
			this->открытьToolStripButton->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->открытьToolStripButton->Name = L"открытьToolStripButton";
			this->открытьToolStripButton->Size = System::Drawing::Size(29, 24);
			this->открытьToolStripButton->Text = L"&Открыть";
			// 
			// сохранитьToolStripButton
			// 
			this->сохранитьToolStripButton->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->сохранитьToolStripButton->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"сохранитьToolStripButton.Image")));
			this->сохранитьToolStripButton->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->сохранитьToolStripButton->Name = L"сохранитьToolStripButton";
			this->сохранитьToolStripButton->Size = System::Drawing::Size(29, 24);
			this->сохранитьToolStripButton->Text = L"&Сохранить";
			// 
			// печатьToolStripButton
			// 
			this->печатьToolStripButton->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->печатьToolStripButton->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"печатьToolStripButton.Image")));
			this->печатьToolStripButton->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->печатьToolStripButton->Name = L"печатьToolStripButton";
			this->печатьToolStripButton->Size = System::Drawing::Size(29, 24);
			this->печатьToolStripButton->Text = L"&Печать";
			// 
			// toolStripSeparator
			// 
			this->toolStripSeparator->Name = L"toolStripSeparator";
			this->toolStripSeparator->Size = System::Drawing::Size(6, 27);
			// 
			// вырезатьToolStripButton
			// 
			this->вырезатьToolStripButton->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->вырезатьToolStripButton->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"вырезатьToolStripButton.Image")));
			this->вырезатьToolStripButton->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->вырезатьToolStripButton->Name = L"вырезатьToolStripButton";
			this->вырезатьToolStripButton->Size = System::Drawing::Size(29, 24);
			this->вырезатьToolStripButton->Text = L"В&ырезать";
			// 
			// копироватьToolStripButton
			// 
			this->копироватьToolStripButton->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->копироватьToolStripButton->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"копироватьToolStripButton.Image")));
			this->копироватьToolStripButton->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->копироватьToolStripButton->Name = L"копироватьToolStripButton";
			this->копироватьToolStripButton->Size = System::Drawing::Size(29, 24);
			this->копироватьToolStripButton->Text = L"&Копировать";
			// 
			// вставкаToolStripButton
			// 
			this->вставкаToolStripButton->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->вставкаToolStripButton->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"вставкаToolStripButton.Image")));
			this->вставкаToolStripButton->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->вставкаToolStripButton->Name = L"вставкаToolStripButton";
			this->вставкаToolStripButton->Size = System::Drawing::Size(29, 24);
			this->вставкаToolStripButton->Text = L"Вст&авка";
			// 
			// toolStripSeparator1
			// 
			this->toolStripSeparator1->Name = L"toolStripSeparator1";
			this->toolStripSeparator1->Size = System::Drawing::Size(6, 27);
			// 
			// справкаToolStripButton
			// 
			this->справкаToolStripButton->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->справкаToolStripButton->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"справкаToolStripButton.Image")));
			this->справкаToolStripButton->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->справкаToolStripButton->Name = L"справкаToolStripButton";
			this->справкаToolStripButton->Size = System::Drawing::Size(29, 24);
			this->справкаToolStripButton->Text = L"Спр&авка";
			// 
			// toolStripButton1
			// 
			this->toolStripButton1->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"toolStripButton1.Image")));
			this->toolStripButton1->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->toolStripButton1->Name = L"toolStripButton1";
			this->toolStripButton1->Size = System::Drawing::Size(69, 24);
			this->toolStripButton1->Text = L"Open";
			this->toolStripButton1->Click += gcnew System::EventHandler(this, &MyForm::toolStripButton1_Click);
			// 
			// toolStripButton2
			// 
			this->toolStripButton2->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"toolStripButton2.Image")));
			this->toolStripButton2->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->toolStripButton2->Name = L"toolStripButton2";
			this->toolStripButton2->Size = System::Drawing::Size(64, 24);
			this->toolStripButton2->Text = L"Save";
			this->toolStripButton2->Click += gcnew System::EventHandler(this, &MyForm::toolStripButton2_Click);
			// 
			// toolStripButton3
			// 
			this->toolStripButton3->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"toolStripButton3.Image")));
			this->toolStripButton3->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->toolStripButton3->Name = L"toolStripButton3";
			this->toolStripButton3->Size = System::Drawing::Size(63, 24);
			this->toolStripButton3->Text = L"Print";
			this->toolStripButton3->Click += gcnew System::EventHandler(this, &MyForm::toolStripButton3_Click);
			// 
			// toolStripButton4
			// 
			this->toolStripButton4->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"toolStripButton4.Image")));
			this->toolStripButton4->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->toolStripButton4->Name = L"toolStripButton4";
			this->toolStripButton4->Size = System::Drawing::Size(87, 24);
			this->toolStripButton4->Text = L"LoadInp";
			this->toolStripButton4->Click += gcnew System::EventHandler(this, &MyForm::toolStripButton4_Click);
			// 
			// toolStripButton5
			// 
			this->toolStripButton5->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"toolStripButton5.Image")));
			this->toolStripButton5->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->toolStripButton5->Name = L"toolStripButton5";
			this->toolStripButton5->Size = System::Drawing::Size(90, 24);
			this->toolStripButton5->Text = L"LoadOut";
			this->toolStripButton5->Click += gcnew System::EventHandler(this, &MyForm::toolStripButton5_Click);
			// 
			// toolStripButton6
			// 
			this->toolStripButton6->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"toolStripButton6.Image")));
			this->toolStripButton6->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->toolStripButton6->Name = L"toolStripButton6";
			this->toolStripButton6->Size = System::Drawing::Size(55, 24);
			this->toolStripButton6->Text = L"Cut";
			this->toolStripButton6->Click += gcnew System::EventHandler(this, &MyForm::toolStripButton6_Click);
			// 
			// toolStripButton7
			// 
			this->toolStripButton7->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"toolStripButton7.Image")));
			this->toolStripButton7->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->toolStripButton7->Name = L"toolStripButton7";
			this->toolStripButton7->Size = System::Drawing::Size(67, 24);
			this->toolStripButton7->Text = L"Copy";
			this->toolStripButton7->Click += gcnew System::EventHandler(this, &MyForm::toolStripButton7_Click);
			// 
			// toolStripButton8
			// 
			this->toolStripButton8->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"toolStripButton8.Image")));
			this->toolStripButton8->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->toolStripButton8->Name = L"toolStripButton8";
			this->toolStripButton8->Size = System::Drawing::Size(67, 24);
			this->toolStripButton8->Text = L"Paste";
			this->toolStripButton8->Click += gcnew System::EventHandler(this, &MyForm::toolStripButton8_Click);
			// 
			// toolStripButton9
			// 
			this->toolStripButton9->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"toolStripButton9.Image")));
			this->toolStripButton9->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->toolStripButton9->Name = L"toolStripButton9";
			this->toolStripButton9->Size = System::Drawing::Size(61, 24);
			this->toolStripButton9->Text = L"Calc";
			this->toolStripButton9->Click += gcnew System::EventHandler(this, &MyForm::toolStripButton9_Click);
			// 
			// richTextBox1
			// 
			this->richTextBox1->Location = System::Drawing::Point(370, 30);
			this->richTextBox1->Name = L"richTextBox1";
			this->richTextBox1->Size = System::Drawing::Size(312, 386);
			this->richTextBox1->TabIndex = 1;
			this->richTextBox1->Text = L"";
			// 
			// treeView1
			// 
			this->treeView1->Location = System::Drawing::Point(29, 30);
			this->treeView1->Name = L"treeView1";
			treeNode1->Name = L"Узел4";
			treeNode1->Tag = L"MLE_Normal";
			treeNode1->Text = L"Нормальный закон распределения";
			treeNode2->Name = L"Узел5";
			treeNode2->Tag = L"MLE_Weibull";
			treeNode2->Text = L"Распределение Вейбулла";
			treeNode3->Name = L"Узел2";
			treeNode3->Text = L"Метод максимального правдоподобия";
			treeNode4->Name = L"Узел0";
			treeNode4->Tag = L"MLS_Normal";
			treeNode4->Text = L"Нормальный закон";
			treeNode5->Name = L"Узел1";
			treeNode5->Tag = L"MLS_Weibull";
			treeNode5->Text = L"Вейбулл";
			treeNode6->Name = L"Узел3";
			treeNode6->Text = L"Метод наименьших квадратов";
			treeNode7->Name = L"Узел2";
			treeNode7->Tag = L"MLS_Regress";
			treeNode7->Text = L"Регрессивный анализ";
			treeNode8->Name = L"Узел1";
			treeNode8->Text = L"Статистическое оценивание";
			treeNode9->Name = L"Узел0";
			treeNode9->Text = L"Обработка данных на ЭВМ";
			this->treeView1->Nodes->AddRange(gcnew cli::array< System::Windows::Forms::TreeNode^  >(1) { treeNode9 });
			this->treeView1->Size = System::Drawing::Size(335, 386);
			this->treeView1->TabIndex = 2;
			// 
			// chart1
			// 
			chartArea1->Name = L"ChartArea1";
			this->chart1->ChartAreas->Add(chartArea1);
			legend1->Name = L"Legend1";
			this->chart1->Legends->Add(legend1);
			this->chart1->Location = System::Drawing::Point(688, 30);
			this->chart1->Name = L"chart1";
			series1->ChartArea = L"ChartArea1";
			series1->Legend = L"Legend1";
			series1->Name = L"Series1";
			this->chart1->Series->Add(series1);
			this->chart1->Size = System::Drawing::Size(807, 692);
			this->chart1->TabIndex = 3;
			this->chart1->Text = L"chart1";
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(8, 16);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1478, 732);
			this->Controls->Add(this->chart1);
			this->Controls->Add(this->treeView1);
			this->Controls->Add(this->richTextBox1);
			this->Controls->Add(this->toolStrip1);
			this->Name = L"MyForm";
			this->Text = L"MyForm";
			this->Load += gcnew System::EventHandler(this, &MyForm::MyForm_Load);
			this->toolStrip1->ResumeLayout(false);
			this->toolStrip1->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void MyForm_Load(System::Object^ sender, System::EventArgs^ e) {
	}
		   void MLE_Normal(string ff) {

			   int i, j, n, k, nc, kx, icount, numres, ifault, kp;
			   const static int nx = 2;
			   int* r, * m, * gtype;
			   string* gtext, * gcolor;
			   string s1;
			   double* x, * y, * xsimpl, * xx, ** v, * fcum, * ycum, * xp, * p;
			   double cp, cko, q, eps, zp, tlow, tup;

			   String^ fff;
			   fff = gcnew String(ff.c_str());

			   ifstream inp("Inp/" + ff + ".inp");
			   ofstream out("Out/" + ff + ".out");
			   inp >> s1; //text
			   inp >> n; //sample size
			   inp >> s1; //text
			   x = new double[n];
			   xsimpl = new double[2];
			   r = new int[n];
			   for (i = 0; i < n; i++) inp >> x[i]; //vector x
			   inp >> s1; //text
			   for (i = 0; i < n; i++) inp >> r[i]; //vector censorized value
			   inp >> s1;
			   inp >> kp;
			   p = new double[kp];
			   inp >> s1;
			   for (i = 0; i < kp; i++) inp >> p[i];
			   inp.close();

			   v = new double* [nx];

			   for (i = 0; i < nx; i++) v[i] = new double[nx];
			   for (i = 0; i < nx; i++) {
				   for (j = 0; j < nx; j++) v[i][j] = 0;
			   }

			   cp = 0; cko = 0; k = 0;
			   for (i = 0; i < n; i++) {
				   k += (1 - r[i]); // количество наблюдений
				   cp += (1 - r[i]) * x[i];
				   cko += (1 - r[i]) * x[i] * x[i];
			   }
			   cp /= k; //выборочное среднее по наблюдениям
			   cko = sqrt((cko - cp * cp * k) / (k - 1)); //выборочное ско по наблюдениям
			   if (k == n) {
				   sm.cp = cp; sm.cko = cko; sm.n = n;
				   q = 0; icount = 0; numres = 0; ifault = 0;
				   xsimpl[0] = cp; xsimpl[1] = cko;
			   }
			   else {
				   sm.cp = cp; sm.cko = cko; sm.n = n;
				   xsimpl[0] = cp; xsimpl[1] = cko; q = 0;
				   sm.x = x; sm.r = r;
				   array<double, nx>start = { sm.cp,sm.cko };
				   array<double, nx>step = { 0.1,0.05 };
				   eps = 1.0e-7;

				   /*
				   lim = 1000; step = 0.01; eps = 0.001;
				   sm.cp = cp; sm.cko = cko; sm.n = n;
				   xsimpl[0] = cp; xsimpl[1] = cko; nx = 2; q = 0;
				   sm.x = x; sm.r = r;
				   simpl(xsimpl, step, eps, lim, ier, nx, q, fun1);
				   */

				   nelder_mead_result<double, nx>result = nelder_mead<double, nx>(Normal_function_to_minimize, start, eps, step);
				   xsimpl[0] = result.xmin[0];
				   xsimpl[1] = result.xmin[1];
				   q = result.ynewlo;
				   icount = result.icount;
				   numres = result.numres;
				   ifault = result.ifault;
				   //simpl(xsimpl, step, eps, lim, ier, nx, q, fun1);



				   CovMatrixMleN(n, x, r, xsimpl[0], xsimpl[1], v);
			   }

			   nc = 4; kx = kp;
			   gtext = new string[nc];
			   gcolor = new string[nc];
			   gtype = new int[nc];
			   gtext[0] = "Distribution"; gtext[1] = "Low"; gtext[2] = "Up"; gtext[3] = "Data";
			   gcolor[0] = "Red"; gcolor[1] = "Green"; gcolor[2] = "Brown"; gcolor[3] = "Black";
			   gtype[0] = 0; gtype[1] = 0; gtype[2] = 0; gtype[3] = 1;

			   m = new int[nc];
			   m[0] = kx; m[1] = kx; m[2] = kx; m[3] = k;
			   y = new double[500];
			   xx = new double[500];
			   xp = new double[2];
			   for (i = 0; i < kx; i++) {
				   zp = invnormaldistribution(p[i]);
				   xx[i] = xsimpl[0] + xsimpl[1] * zp;
				   y[i] = zp + 5.0;
				   y[i + kx] = zp + 5;
				   y[i + 2 * kx] = zp + 5;
				   if (k == n) {
					   lmtexact(0.05, n, p[i], tlow);
					   lmtexact(0.95, n, p[i], tup);
				   }
				   else {
					   lmtaprn(0.95, n, v[0][0], v[1][1], v[0][1], zp * sqrt(n), xp);
					   tlow = xp[0]; tup = xp[1];
				   }
				   xx[i + kx] = xsimpl[0] + (tlow * xsimpl[1]) / sqrt(n);
				   xx[i + 2 * kx] = xsimpl[0] + (tup * xsimpl[1]) / sqrt(n);
			   }
			   fcum = new double[n];
			   ycum = new double[n];
			   cum(n, x, r, k, fcum, ycum);
			   for (i = 0; i < k; i++) {
				   xx[i + 3 * kx] = ycum[i];
				   y[i + 3 * kx] = 5.0 + invnormaldistribution(fcum[i]);
			   }

			   Graph(nc, m, xx, y, gtype, gtext, gcolor);


			   out << "Method:" << ff << "\n";
			   out << "n=" << n << "\n";
			   out << "X" << "\n";
			   for (i = 0; i < n; i++) out << x[i] << " , ";
			   out << "\n";
			   out << "R" << "\n";
			   for (i = 0; i < n; i++) out << r[i] << " , ";
			   out << "\n";
			   out << "cp*=" << sm.cp << "\n";
			   out << "cko*=" << sm.cko << "\n";
			   out << "Q=" << q << "\n";
			   out << "icount=" << icount << endl;
			   out << "numres=" << numres << endl;
			   out << "ifault=" << ifault << endl;
			   out << "x[0]=" << xsimpl[0] << "\n";
			   out << "x[1]=" << xsimpl[1] << "\n";

			   out << "P" << "\n";
			   for (i = 0; i < kx; i++) out << p[i] << " ; ";
			   out << "\n";
			   out << "Xplow" << "\n";
			   for (i = 0; i < kx; i++) out << xx[i + kx] << " ; ";
			   out << "\n";
			   out << "Xp" << "\n";
			   for (i = 0; i < kx; i++) out << xx[i] << " ; ";
			   out << "\n";
			   out << "Xpup" << "\n";
			   for (i = 0; i < kx; i++) out << xx[i + 2 * kx] << " ; ";
			   out << "\n";

			   out << "v11=" << v[0][0] << "\n";
			   out << "v12=" << v[0][1] << "\n";
			   out << "v21=" << v[1][0] << "\n";
			   out << "v22=" << v[1][1] << "\n";
			   out.close();
			   richTextBox1->Clear();
			   richTextBox1->Text = File::ReadAllText("Out/" + fff + ".out");
			   delete[] xsimpl, xx, x, y, fcum, ycum, v, r, m, gtext, gcolor, gtype, xp, p;
		   }
	private: System::Void toolStripButton1_Click(System::Object^ sender, System::EventArgs^ e) {
		wchar_t path[MAX_PATH];
		GetCurrentDirectory(sizeof(path), path);

		OpenFileDialog^ openFileDialog1 = gcnew OpenFileDialog;
		openFileDialog1->InitialDirectory = Convert::ToString(path);
		openFileDialog1->Filter = "All Files|*.*";
		openFileDialog1->FilterIndex = 2;
		openFileDialog1->RestoreDirectory = true;
		if (openFileDialog1->ShowDialog() == System::Windows::Forms::DialogResult::OK) {
			richTextBox1->Text = File::ReadAllText(openFileDialog1->FileName);
		}
	}
	private: System::Void toolStripButton2_Click(System::Object^ sender, System::EventArgs^ e) {
		wchar_t path[MAX_PATH];

		GetCurrentDirectory(sizeof(path), path);
		SaveFileDialog^ saveFileDialog = gcnew SaveFileDialog();
		saveFileDialog->InitialDirectory = Convert::ToString(path);
		saveFileDialog->Filter = "All Files|*.*";
		if (saveFileDialog->ShowDialog() == System::Windows::Forms::DialogResult::OK) {
			File::WriteAllText(saveFileDialog->FileName, richTextBox1->Text);
		}
	}
	private: System::Void toolStripButton3_Click(System::Object^ sender, System::EventArgs^ e) {
		wchar_t path[MAX_PATH];
		GetCurrentDirectory(sizeof(path), path);
		PrintDialog^ printDialog = gcnew PrintDialog();
		printDialog->ShowDialog();
		if (printDialog->ShowDialog() == System::Windows::Forms::DialogResult::OK) {

		}
	}
	private: System::Void toolStripButton4_Click(System::Object^ sender, System::EventArgs^ e) {
		String^ ff;

		ff = Convert::ToString(treeView1->SelectedNode->Tag);
		if (ff == "") return;
		ff = ff + ".inp";
		richTextBox1->Text = File::ReadAllText(ff);
	}
	private: System::Void toolStripButton5_Click(System::Object^ sender, System::EventArgs^ e) {
		String^ ff;
		ff = Convert::ToString(treeView1->SelectedNode->Tag);
		if (ff == "") return;
		ff = ff + ".out";
		richTextBox1->Text = File::ReadAllText(ff);
	}
	private: System::Void toolStripButton6_Click(System::Object^ sender, System::EventArgs^ e) {
		if (richTextBox1->SelectionLength > 0) richTextBox1->Cut();
	}
	private: System::Void toolStripButton7_Click(System::Object^ sender, System::EventArgs^ e) {
		if (richTextBox1->SelectionLength > 0) richTextBox1->Copy();
	}
	private: System::Void toolStripButton8_Click(System::Object^ sender, System::EventArgs^ e) {
		if (Clipboard::GetDataObject()->GetDataPresent(DataFormats::Rtf)) { //Есть Rtf в буфере
			if (richTextBox1->SelectionLength > 0) { //И что-то выделено,
			//спросим, как вставлять - поверх выделенного или в конец?
				if (MessageBox::Show(L"Вставить поверх выделения?", L"Сообщение",
					MessageBoxButtons::YesNo) == System::Windows::Forms::DialogResult::No)
					richTextBox1->SelectionStart = richTextBox1->Text->Length;
			}
			richTextBox1->Paste();
		}
	}
	private: System::Void toolStripButton9_Click(System::Object^ sender, System::EventArgs^ e) {
		string ff;
		String^ fff;
		fff = Convert::ToString(treeView1->SelectedNode->Tag);
		if (fff == "") return;

		ff = marshal_as<string>(fff);
		if (fff == "MLE_Normal") {
			MLE_Normal(ff);
			return;
		}
		if (fff == "MLE_Weibull") {
			MLE_Weibull(ff);
			return;
		}
		if (ff == "MLS_Normal" || ff == "MLS_Weibull") {
			MLS(ff);
			return;
		}
		if (ff == "MLS_Regress") {
			MLS_Regress(ff);
			return;

		}
	}
		   //#######################################################################
		   void MLE_Weibull(string ff) {
			   int i, j, n, k, nc, kx, kp;
			   const static int nx = 1;
			   int* r, * m, * gtype;
			   string* gtext, * gcolor;
			   string s1;
			   double* x, * y, * xx, * xsimpl, ** v, * fcum, * ycum, * xp, * p;
			   double cp, cko, eps, q, zp, s;

			   String^ fff;

			   fff = gcnew String(ff.c_str());


			   ifstream inp("Inp/" + ff + ".inp");
			   inp >> s1;
			   inp >> n;
			   inp >> s1;
			   x = new double[n];
			   r = new int[n];
			   for (i = 0; i < n; i++) inp >> x[i];
			   inp >> s1;
			   for (i = 0; i < n; i++) inp >> r[i];
			   inp >> s1;
			   inp >> kp;
			   p = new double[kp];
			   inp >> s1;
			   for (i = 0; i < kp; i++) inp >> p[i];
			   inp.close();

			   v = new double* [2];
			   for (i = 0; i < 2; i++) v[i] = new double[2];
			   for (i = 0; i < 2; i++) {
				   for (j = 0; j < 2; j++) v[i][j] = 0;
			   }

			   y = new double[500];
			   xx = new double[500];
			   xsimpl = new double[2];
			   xp = new double[2];

			   q = 0;
			   cp = 2; cko = 0; sm.cp = cp; sm.cko = cko; sm.n = n;
			   sm.x = x; sm.r = r;
			   array<double, nx>start = { 2. };
			   array<double, nx>step = { 0.1 };
			   eps = 1.0e-12;
			   /*
			   lim = 10000; step = 0.1; eps = 0.00001; q = 0;
			   nx = 1; xsimpl[0] = 2.; xsimpl[1] = 0;
			   cp = 2; cko = 0; sm.cp = cp; sm.cko = cko; sm.n = n;
			   sm.x = x; sm.r = r;
			   simpl(xsimpl, step, eps, lim, ier, nx, q, fun2);
			   s = 0.; k = 0;
			   for (i = 0; i < n; i++) { k += (1 - r[i]); s += pow(x[i], xsimpl[0]); }
			   xsimpl[1] = pow(s / k, 1 / xsimpl[0]);
			   */
			   nelder_mead_result<double, nx>result = nelder_mead<double, nx>(Weibull_function_to_minimize, start, eps, step);
			   xsimpl[0] = result.xmin[0];
			   q = result.ynewlo;

			   s = 0.; k = 0;
			   for (i = 0; i < n; i++) { k += (1 - r[i]); s += pow(x[i], xsimpl[0]); }
			   xsimpl[1] = pow(s / k, 1 / xsimpl[0]);

			   CovMatrixMleW(n, x, r, xsimpl[1], xsimpl[0], v);

			   nc = 4; kx = kp;
			   gtext = new string[nc];
			   gcolor = new string[nc];
			   gtype = new int[nc];
			   gtext[0] = "Distribution"; gtext[1] = "Low"; gtext[2] = "Up"; gtext[3] = "Data";
			   gcolor[0] = "Black"; gcolor[1] = "Green"; gcolor[2] = "Brown"; gcolor[3] = "Red";
			   gtype[0] = 0; gtype[1] = 0; gtype[2] = 0; gtype[3] = 1;
			   m = new int[nc];
			   m[0] = kx; m[1] = kx; m[2] = kx; m[3] = k;

			   xsimpl[0] = 1 / xsimpl[0]; xsimpl[1] = log(xsimpl[1]);
			   for (i = 0; i < kx; i++) {
				   zp = log(log(1 / (1 - p[i])));
				   xx[i] = xsimpl[1] + xsimpl[0] * zp; y[i] = zp;
				   lmtaprn(0.95, n, v[0][0], v[1][1], v[0][1], zp * sqrt(n), xp);
				   xx[i + kx] = xsimpl[1] + (xp[0] * xsimpl[0]) / sqrt(n); y[i + kx] = zp;
				   xx[i + 2 * kx] = xsimpl[1] + (xp[1] * xsimpl[0]) / sqrt(n); y[i + 2 * kx] = zp;
			   }
			   fcum = new double[n];
			   ycum = new double[n];
			   cum(n, x, r, k, fcum, ycum);
			   for (i = 0; i < k; i++) {
				   xx[i + 3 * kx] = log(ycum[i]);
				   y[i + 3 * kx] = log(log(1 / (1 - fcum[i])));
			   }
			   Graph(nc, m, xx, y, gtype, gtext, gcolor);

			   ofstream out("Out/" + ff + ".out");
			   out << "Method:" << ff << "\n";
			   out << "n=" << n << "\n";
			   out << "X" << "\n";
			   for (i = 0; i < n; i++) out << x[i] << " , ";
			   out << "\n";
			   out << "R" << "\n";
			   for (i = 0; i < n; i++) out << r[i] << " , ";
			   out << "\n";
			   out << "cp*=" << sm.cp << "\n";
			   out << "cko*=" << sm.cko << "\n";
			   out << "Q=" << q << "\n";
			   out << "icount=" << result.icount << endl;
			   out << "numres=" << result.numres << endl;
			   out << "ifault=" << result.ifault << endl;
			   out << "x[0]=" << xsimpl[0] << "\n";
			   out << "x[1]=" << xsimpl[1] << "\n";

			   out << "P" << "\n";
			   for (i = 0; i < kx; i++) out << p[i] << " ; ";
			   out << "\n";
			   out << "Xplow" << "\n";
			   for (i = 0; i < kx; i++) out << xx[i + kx] << " ; ";
			   out << "\n";
			   out << "Xp" << "\n";
			   for (i = 0; i < kx; i++) out << xx[i] << " ; ";
			   out << "\n";
			   out << "Xpup" << "\n";
			   for (i = 0; i < kx; i++) out << xx[i + 2 * kx] << " ; ";
			   out << "\n";

			   out << "v11=" << v[0][0] << "\n";
			   out << "v12=" << v[0][1] << "\n";
			   out << "v21=" << v[1][0] << "\n";
			   out << "v22=" << v[1][1] << "\n";
			   out.close();
			   richTextBox1->Clear();
			   richTextBox1->Text = File::ReadAllText("Out/" + fff + ".out");
			   delete[] xsimpl, xx, x, y, fcum, ycum, v, r, m, gtext, gcolor, gtype, xp, p;
		   }
		   //##########################################################################
		   // ##########Chart####################################################
			   void Graph(int nc, int* m, double* x, double* y, int* gtype, string * gtext, string * gcolor) {
			   int i, j, kx;
			   double z, f;
			   using namespace System::Collections::Generic;
			   using namespace System::Drawing::Drawing2D;
			   using namespace System::Windows::Forms::DataVisualization::Charting;

			   chart1->Series->Clear();
			   kx = 0;
			   for (i = 0; i < nc; i++) {
				   chart1->Series->Add("mySeries" + Convert::ToString(i));
				   if (gtype[i] == 1) {
					   chart1->Series[i]->ChartType = SeriesChartType::Point;
					   chart1->Series[i]->Color = System::Drawing::Color::FromName(gcnew System::String(gcolor[i].c_str()));
					   chart1->Series[i]->MarkerStyle = MarkerStyle::Star6; // Diamond;Square;Triangle;Cross;Circle;
					   chart1->Series[i]->MarkerColor = System::Drawing::Color::FromName(gcnew System::String(gcolor[i].c_str()));
					   chart1->Series[i]->MarkerSize = 10;
				   }
				   if (gtype[i] == 0) {
					   chart1->Series[i]->ChartType = SeriesChartType::Line;
					   chart1->Series[i]->Color = System::Drawing::Color::FromName(gcnew System::String(gcolor[i].c_str()));
					   chart1->Series[i]->BorderWidth = 3;
				   }
				   chart1->Series[i]->LegendText = gcnew System::String(gtext[i].c_str());

				   for (j = 0; j < m[i]; j++) {
					   z = x[j + kx];
					   z = floor(z * 100) / 100;
					   f = y[j + kx];
					   chart1->Series[i]->Points->AddXY(z, f);
				   }
				   kx += m[i];
			   }
		   }
			   void MLS(string ff) {

				   int i, j, n, k, nc, kx, ku, kp;
				   int* m, * gtype;
				   string* gtext, * gcolor;
				   string s1;
				   double* xz, ** y, ** x, * yy, * xp, * xx, ** v, * yr, ** b, ** db, * p;
				   double cp, cko, zp, tlow, tup, vrs, er;

				   String^ fff;

				   if (ff == "MLS_Normal") ku = 0;
				   if (ff == "MLS_Weibull") ku = 1;

				   fff = gcnew String(ff.c_str());
				   ifstream inp("Inp/" + ff + ".inp");
				   inp >> s1;
				   inp >> n;
				   inp >> s1;
				   xz = new double[n];
				   for (i = 0; i < n; i++) inp >> xz[i];
				   inp >> s1;
				   inp >> kp;
				   p = new double[kp];
				   inp >> s1;
				   for (i = 0; i < kp; i++) inp >> p[i];
				   inp.close();


				   k = 2;
				   yr = new double[n];
				   x = new double* [n];
				   y = new double* [n];
				   v = new double* [n];
				   db = new double* [k];
				   b = new double* [k];

				   for (i = 0; i < n; i++) { x[i] = new double[n]; y[i] = new double[n]; v[i] = new double[n]; }
				   for (i = 0; i < k; i++) { db[i] = new double[k]; b[i] = new double[k]; }
				   for (i = 0; i < k; i++) {
					   for (j = 0; j < k; j++) { b[i][j] = 0; db[i][j] = 0; }
				   }

				   for (i = 0; i < n; i++) {
					   if (ku == 1) xz[i] = log(xz[i]);
					   y[i][0] = xz[i];
				   }

				   standart(n, xz, cp, cko);

				   for (i = 0; i < n; i++) {
					   for (j = i; j < n; j++) {
						   if (ku == 0) ordern(n, (i + 1.) / (n + 1.), (j + 1) / (n + 1.), er, vrs);
						   if (ku == 1) orderw(n, (i + 1.) / (n + 1.), (j + 1) / (n + 1.), er, vrs);
						   v[j][i] = vrs; v[i][j] = vrs;
					   }
					   x[i][0] = 1; x[i][1] = er;
				   }


				   MleastSquare_weight(n, k, x, y, v, db, b, yr);

				   xx = new double[500];
				   yy = new double[500];
				   xp = new double[2];

				   nc = 4;
				   kx = kp;

				   gtext = new string[nc];
				   gcolor = new string[nc];
				   gtype = new int[nc];

				   gtext[0] = "Distribution"; gtext[1] = "Low"; gtext[2] = "Up"; gtext[3] = "Data";
				   gcolor[0] = "Black"; gcolor[1] = "Green"; gcolor[2] = "Brown"; gcolor[3] = "Red";
				   gtype[0] = 0; gtype[1] = 0; gtype[2] = 0; gtype[3] = 1;

				   m = new int[nc];
				   m[0] = kx; m[1] = kx; m[2] = kx; m[3] = n;
				   for (i = 0; i < kx; i++) {
					   if (ku == 0) {
						   zp = invnormaldistribution(p[i]);
						   lmtexact(0.05, n, p[i], tlow);
						   lmtexact(0.95, n, p[i], tup);
					   }
					   if (ku == 1) {
						   zp = log(log(1 / (1 - p[i])));
						   lmtaprn(0.95, n, db[0][0], db[1][1], db[0][1], zp * sqrt(n), xp);
						   tlow = xp[0]; tup = xp[1];
					   }
					   xx[i] = b[0][0] + b[1][0] * zp;
					   yy[i] = zp + 5.0;
					   xx[i + kx] = b[0][0] + (tlow * b[1][0]) / sqrt(n); yy[i + kx] = zp + 5;
					   xx[i + 2 * kx] = b[0][0] + (tup * b[1][0]) / sqrt(n); yy[i + 2 * kx] = zp + 5;
				   }
				   for (i = 0; i < n; i++) {
					   if (ku == 0) zp = invnormaldistribution((i + 1.) / (n + 1.));
					   if (ku == 1) zp = log(log(1 / (1 - (i + 1.) / (n + 1.))));
					   xx[i + 3 * kx] = xz[i];
					   yy[i + 3 * kx] = 5.0 + zp;
				   }
				   Graph(nc, m, xx, yy, gtype, gtext, gcolor);

				   ofstream out("Out/" + ff + ".out");
				   out << "Method:" << ff << "\n";
				   out << "n=" << n << "\n";
				   out << "X" << "\n";
				   for (i = 0; i < n; i++) out << xz[i] << " , ";
				   out << "\n";
				   out << "Yr" << "\n";
				   for (i = 0; i < n; i++) out << yr[i] << " , ";
				   out << "\n";
				   out << "cp*=" << cp << "\n";
				   out << "cko*=" << cko << "\n";
				   out << "cp=" << b[0][0] << "\n";
				   out << "cko=" << b[1][0] << "\n";

				   out << "P" << "\n";
				   for (i = 0; i < kx; i++) out << p[i] << " ; ";
				   out << "\n";
				   out << "Xplow" << "\n";
				   for (i = 0; i < kx; i++) out << xx[i + kx] << " ; ";
				   out << "\n";
				   out << "Xp" << "\n";
				   for (i = 0; i < kx; i++) out << xx[i] << " ; ";
				   out << "\n";
				   out << "Xpup" << "\n";
				   for (i = 0; i < kx; i++) out << xx[i + 2 * kx] << " ; ";
				   out << "\n";

				   out << "t11=" << db[0][0] << "\n";
				   out << "t22=" << db[0][1] << "\n";
				   out << "t12=" << db[1][1] << "\n";
				   out.close();
				   richTextBox1->Clear();
				   richTextBox1->Text = File::ReadAllText("Out/" + fff + ".out");
				   delete[] xz, xx, yy, x, y, yr, v, m, gtext, gcolor, gtype, xp, p;
			   }
			   //######################################################################
			   void MLS_Regress(string ff) {
				   string s1;
				   int i, j, k, n, nc;
				   double** x, ** y, ** db, ** b, * yr, * xz, * yz;
				   string* gtext, * gcolor;
				   int* m, * gtype;
				   double s;
				   String^ fff;

				   fff = gcnew String(ff.c_str());
				   ifstream inp("Inp/" + ff + ".inp");
				   ofstream out("Out/" + ff + ".out");
				   inp >> s1;
				   inp >> k;
				   inp >> s1;
				   inp >> n;

				   xz = new double[500];
				   yz = new double[500];
				   x = new double* [n];
				   y = new double* [n];
				   db = new double* [k];
				   b = new double* [k];
				   yr = new double[n];
				   for (i = 0; i < n; i++) { x[i] = new double[n]; y[i] = new double[n]; }
				   for (i = 0; i < k; i++) { db[i] = new double[k]; b[i] = new double[k]; }


				   for (i = 0; i < k; i++) {
					   for (j = 0; j < k; j++) db[i][j] = 0;
				   }
				   for (i = 0; i < k; i++) b[i][0] = 0;


				   inp >> s1;
				   for (i = 0; i < n; i++) {
					   for (j = 0; j < k; j++) inp >> x[i][j];
				   }
				   inp >> s1;
				   for (i = 0; i < n; i++) inp >> y[i][0];
				   inp.close();

				   db = InverseMatrix((MultiplyMatrix(k, n, n, k, TransMatrix(n, k, x), x)), k); //covariance matrix factors (k x k)
				   b = MultiplyMatrix(k, n, n, 1, MultiplyMatrix(k, k, k, n, db, TransMatrix(n, k, x)), y); // coef

				   out << "k=" << k << "\n";
				   out << "n=" << n << "\n";
				   out << "X" << "\n";
				   for (i = 0; i < n; i++) {
					   for (j = 0; j < k; j++) out << x[i][j] << " , ";
					   out << "\n";
				   }
				   out << "Y" << "\n";
				   for (i = 0; i < n; i++) out << y[i][0] << " , ";
				   out << "\n";
				   out << "b" << "\n";
				   for (i = 0; i < k; i++) out << b[i][0] << " , ";
				   out << "\n";
				   out << "D{b}" << "\n";
				   for (i = 0; i < k; i++) {
					   for (j = 0; j < k; j++) out << db[i][j] << " , ";
					   out << "\n";
				   }

				   for (i = 0; i < n; i++) {
					   s = 0;
					   for (j = 0; j < k; j++) s += b[j][0] * x[i][j];
					   yr[i] = s;
				   }
				   out << "Yr" << "\n";
				   for (i = 0; i < n; i++) out << yr[i] << " , ";
				   out.close();

				   richTextBox1->Clear();
				   richTextBox1->Text = File::ReadAllText("Out/" + fff + ".out");

				   nc = 2;
				   m = new int[nc];
				   gtext = new string[nc];
				   gcolor = new string[nc];
				   gtype = new int[nc];

				   gtext[0] = "Data"; gtext[1] = "Regress Curve";
				   gcolor[0] = "Black"; gcolor[1] = "Red";
				   gtype[0] = 1; gtype[1] = 0;
				   nc = 2;
				   m[0] = n;
				   for (i = 0; i < n; i++) {
					   xz[i] = i + 1.0;
					   yz[i] = y[i][0];
				   }
				   m[1] = n;
				   for (i = 0; i < n; i++) {
					   xz[i + n] = i + 1.0;
					   yz[i + n] = yr[i];
				   }
				   Graph(nc, m, xz, yz, gtype, gtext, gcolor);
				   delete[] x, y, db, b, yr, xz, yz, gtext, gcolor, m, gtype;
			   }
};
}
