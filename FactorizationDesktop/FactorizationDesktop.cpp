// FactorizationDesktop.cpp : Defines the entry point for the application.
//
#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#include <windows.h>

// C RunTime Header Files
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <tchar.h>
#include "FactorizationDesktop.h"
#include <thread>
#include <string>
#include <vector>
#include "Header.h"
#define MAX_LOADSTRING 100

// Global Variables:
HINSTANCE hInst;                                // current instance
WCHAR szTitle[MAX_LOADSTRING];                  // The title bar text
WCHAR szWindowClass[MAX_LOADSTRING];            // the main window class name





HWND TextBox;
HWND ProgressBox;
HWND ResultBox;
HWND mainWindow;
HWND button;

// Forward declarations of functions included in this code module:
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: Place code here.

    // Initialize global strings
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_FACTORIZATIONDESKTOP, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);

    // Perform application initialization:
    if (!InitInstance (hInstance, nCmdShow))
    {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_FACTORIZATIONDESKTOP));

    MSG msg;

    // Main message loop:
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    return (int) msg.wParam;
}



//
//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_FACTORIZATIONDESKTOP));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW+1);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_FACTORIZATIONDESKTOP);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

//
//   FUNCTION: InitInstance(HINSTANCE, int)
//
//   PURPOSE: Saves instance handle and creates main window
//
//   COMMENTS:
//
//        In this function, we save the instance handle in a global variable and
//        create and display the main program window.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   hInst = hInstance; // Store instance handle in our global variable
   
   HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
      CW_USEDEFAULT, 0, 600, 600, nullptr, nullptr, hInstance, nullptr);
   mainWindow = hWnd;
   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, nCmdShow);
   UpdateWindow(hWnd);

   return TRUE;
}

std::wstring s2ws(const std::string& s)
{
	int len;
	int slength = (int)s.length() + 1;
	len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0);
	wchar_t* buf = new wchar_t[len];
	MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
	std::wstring r(buf);
	delete[] buf;
	return r;
}

void sendProgress(std::string progress) {
	std::wstring stemp = s2ws(progress);
	LPCWSTR resultString = stemp.c_str();
	SendMessage(mainWindow, WM_COMMAND, IDM_PROGRESS, (LPARAM)resultString);
}

void threadFunction(HWND hWnd, wchar_t str[])
{
	wchar_t* t = &str[0];
	std::wstring ws(t);
	std::string s(ws.begin(), ws.end());
	if (s.size() == 0) {
		SendMessage(hWnd, WM_COMMAND, IDM_RESULT, (LPARAM)(L"Empty number"));
		EnableWindow(button, true);
		return;
	}
	if (s.size() == 1 && s[0]=='1') {
		SendMessage(hWnd, WM_COMMAND, IDM_RESULT, (LPARAM)(L"1"));
		EnableWindow(button, true);
		return;
	}
	std::vector<std::string> result = fact::fatorize(s);
	uint64_t size = result.size();
	for (uint64_t i = 0; i < size; ++i) {
		std::wstring stemp = s2ws(result[i]);
		LPCWSTR resultString = stemp.c_str();
		SendMessage(hWnd, WM_COMMAND, IDM_RESULT, (LPARAM)resultString);
	}
	EnableWindow(button, true);
}

void AppendText(HWND hEditWnd, LPCTSTR Text)
{
	int idx = GetWindowTextLength(hEditWnd);
	SendMessage(hEditWnd, EM_SETSEL, (WPARAM)idx, (LPARAM)idx);
	SendMessage(hEditWnd, EM_REPLACESEL, 0, (LPARAM)Text);
}
void AppendLine(HWND hEditWnd, LPARAM Text)
{
	LPCWSTR nl = L"\r\n";
	int idx = GetWindowTextLength(hEditWnd);
	SendMessage(hEditWnd, EM_SETSEL, (WPARAM)idx, (LPARAM)idx);
	SendMessage(hEditWnd, EM_REPLACESEL, 0, Text);
	idx = GetWindowTextLength(hEditWnd);
	SendMessage(hEditWnd, EM_SETSEL, (WPARAM)idx, (LPARAM)idx);
	SendMessage(hEditWnd, EM_REPLACESEL, 0, (LPARAM)nl);

}
//
//  FUNCTION: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  PURPOSE:  Processes messages for the main window.
//
//  WM_COMMAND  - process the application menu
//  WM_PAINT    - Paint the main window
//  WM_DESTROY  - post a quit message and return
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    switch (message)
    {
	case WM_CREATE:
		TextBox = CreateWindow(L"EDIT", L"", 
			WS_BORDER | WS_CHILD | WS_VISIBLE | ES_NUMBER, 
			10, 10, 400, 20,
			hWnd, NULL, NULL, NULL);
		ProgressBox = CreateWindow(L"EDIT", L"",
			WS_BORDER | WS_CHILD | WS_VISIBLE | ES_READONLY,
			10, 40, 400, 20,
			hWnd, NULL, NULL, NULL);
		button = CreateWindow(L"BUTTON", L"Factorize", 
			WS_BORDER | WS_CHILD | WS_VISIBLE, 
			420, 10, 70, 20, 
			hWnd, (HMENU)IDM_BUTT0N, NULL, NULL);
		ResultBox = CreateWindow(L"EDIT", L"", 
			WS_BORDER | WS_CHILD | WS_VISIBLE | WS_VSCROLL| ES_LEFT | ES_MULTILINE | ES_READONLY,
			10, 70, 400, 400, 
			hWnd, NULL, NULL, NULL);
		break;
    case WM_COMMAND:
        {
            int wmId = LOWORD(wParam);
            // Parse the menu selections:
            switch (wmId)
            {
			case IDM_BUTT0N: {				
				int gwtstat;
				wchar_t str[255];
				wchar_t* t = &str[0];
				gwtstat = GetWindowText(TextBox, t, 255);
				SetWindowText(ResultBox, L"");
				SetWindowText(ProgressBox, L"");
				std::thread thr(threadFunction, hWnd, str);
				thr.detach();
				EnableWindow(button, false);
				break;
			}
			case IDM_RESULT: {
				AppendLine(ResultBox, lParam);
				break;
			}
			case IDM_PROGRESS: {
				SetWindowText(ProgressBox, (LPCWSTR)lParam);
				break;
			}
            case IDM_ABOUT:
                DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
                break;
            case IDM_EXIT:
                DestroyWindow(hWnd);
                break;
            default:
                return DefWindowProc(hWnd, message, wParam, lParam);
            }
        }
        break;
    case WM_PAINT:
        {
            PAINTSTRUCT ps;
            HDC hdc = BeginPaint(hWnd, &ps);
            // TODO: Add any drawing code that uses hdc here...
            EndPaint(hWnd, &ps);
        }
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}

// Message handler for about box.
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}
