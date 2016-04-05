//
//  StringArray.h
//  CatsReformat
//
//  Created by Jennifer Johnson on 5/2/15.
//  Copyright (c) 2015 Jennifer Johnson. All rights reserved.
//

#ifndef __STRING_ARRAY_H__
#define __STRING_ARRAY_H__

#include "StringBasics.h"

class StringArray
{
protected:
    String ** strings;
    int size, count;
    
public:
    static int alloc;
    
    StringArray(int startsize = 0);
    StringArray(StringArray & original);
    virtual ~StringArray();
    
    // Each line in a file is parsed into a separate array element
    //
    
    void Read(FILE * f);
    void Write(FILE * f);
    void Read(const char * filename);
    void Write(const char * filename);
    
    // Write all strings to the screen
    void Print();
    
    void Grow(int newsize);
    void Clear();
    
    int Length() const { return count; }
    int Dimension(int newcount);
    int CharLength();
    
    String & operator [] (int i) { return *(strings[i]); }
    const String & operator [] (int i) const { return *(strings[i]); }
    
    // These functions divide a string into tokens and append these to the
    // array. Return value is the new array length
    //
    
    int AddColumns(const String & s, char ch = '\t');
    int AddTokens(const String & s, char ch);
    int AddTokens(const String & s, const String & separators = " \t\r\n");
    
    int ReplaceColumns(const String & s, char ch = '\t')
    { Clear(); return AddColumns(s, ch); }
    int ReplaceTokens(const String & s, const String & separators = " \t\r\n")
    { Clear(); return AddTokens(s, separators); }
    
    // These functions add, insert or remove a single array element
    //
    
    int  Add(const String & s);
    void InsertAt(int position, const String & s);
    void Delete(int position);
    
    // These functions manipulate a string as a stack
    //
    
    String & Last() const;
    int      Push(const String & s) { return Add(s); }
    String   Pop();
    
    // Linear search (N/2 comparisons on average) for a single element
    // If searching is required, StringMaps are a better option
    //
    
    int Find(const String & s) const;
    int FastFind(const String & s) const;
    int SlowFind(const String & s) const;
    
    // Alphetically orders strings
    //
    void Sort();
    
    // Trims strings to remove whitespace
    void Trim();
    
    StringArray & operator = (const StringArray & rhs);
    
    bool operator == (const StringArray & rhs);
    bool operator != (const StringArray & rhs)
    { return !(*this == rhs); }
    
private:
    static int ComparisonForSort(const void * a, const void * b);
};

#endif

