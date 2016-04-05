//
//  StringHash.h
//  CatsReformat
//
//  Created by Jennifer Johnson on 5/2/15.
//  Copyright (c) 2015 Jennifer Johnson. All rights reserved.
//

#ifndef __STRINGHASH_H__
#define __STRINGHASH_H__

#include "StringBasics.h"
#include "Constant.h"
#include "Hash.h"

class StringHash
{
protected:
    String        ** strings;
    void          ** objects;
    unsigned int      * keys;
    unsigned int count, size;
    unsigned int        mask;
    
public:
    StringHash(int startsize = 32);
    virtual ~StringHash();
    
    void Grow()    { SetSize(size * 2); }
    void Shrink()  { SetSize(size / 2); }
    
    void SetSize(int newsize);
    
    void Clear();
    
    int  Capacity() const { return size; }
    int  Entries() const  { return count; }
    
    void * Object(int i) const { return objects[i]; }
    void * Object(const String & key) const
    {
        int index = Find(key);
        
        return index >= 0 ? objects[index] : NULL;
    }
    void * Object(const String & key, void * (*create_object)())
    { return objects[Find(key, create_object)]; }
    
    void SetObject(int i, void * object)
    { objects[i] = object; }
    void SetObject(const String & key, void * object)
    { Add(key, object); }
    
    int Add(const String & s, void * object = NULL);
    int Find(const String & s, void * (*create_object)() = NULL);
    int Find(const String & s) const;
    
    StringHash & operator = (const StringHash & rhs);
    
    const String & operator [] (int i) const { return *(strings[i]); }
    String & operator [] (int i) { return *(strings[i]); }
    //      String & String(int i) { return *(strings[i]); }
    
    static void * CreateHash();
    
    void Delete(unsigned int index);
    void Delete(const String & key) { Delete(Find(key)); }
    
    bool SlotInUse(int index) { return strings[index] != NULL; }
    
    void Print();
    void Print(FILE * file);
    void Print(const char * filename);
    
private:
    
    unsigned int Iterate(unsigned int key, const String & string) const
    {
        unsigned int h = key & mask;
        
        while (  strings[h] != NULL &&
               ( keys[h] != key ||
                strings[h]->SlowCompare(string) != 0) )
            h = (h + 1) & mask;
        
        return h;
    }
    
    void Insert(unsigned int where, unsigned int key, const String & string)
    {
        strings[where] = new String;
        *(strings[where]) = string;
        keys[where] = key;
        
        count++;
    }
};

class StringIntHash
{
protected:
    String        ** strings;
    int           * integers;
    unsigned int      * keys;
    unsigned int count, size;
    unsigned int        mask;
    
public:
    StringIntHash(int startsize = 32);
    virtual ~StringIntHash();
    
    void Grow()    { SetSize(size * 2); }
    void Shrink()  { SetSize(size / 2); }
    
    void SetSize(int newsize);
    
    void Clear();
    
    int  Capacity() const { return size; }
    int  Entries() const  { return count; }
    
    int Integer(int i) const { return integers[i]; }
    int Integer(const String & key) const
    {
        int index = Find(key);
        
        return index >= 0 ? integers[index] : -1;
    }
    
    void SetInteger(int i, int value)
    { integers[i] = value; }
    void SetInteger(const String & key, int value)
    { Add(key, value); }
    
    int IncrementCount(const String & key);
    int DecrementCount(const String & key);
    int GetCount(const String & key) const;
    int GetCount(int index) const { return integers[index]; }
    
    int Add(const String & s, int integer);
    int Find(const String & s, int defaultValue);
    int Find(const String & s) const;
    
    StringIntHash & operator = (const StringIntHash & rhs);
    
    const String & operator [] (int i) const { return *(strings[i]); }
    String & operator [] (int i) { return *(strings[i]); }
    //      String & String(int i) { return *(strings[i]); }
    
    void Delete(unsigned int index);
    void Delete(const String & key) { Delete(Find(key)); }
    
    bool SlotInUse(int index) { return strings[index] != NULL; }
    
private:
    
    unsigned int Iterate(unsigned int key, const String & string) const
    {
        unsigned int h = key & mask;
        
        while (  strings[h] != NULL &&
               ( keys[h] != key ||
                strings[h]->SlowCompare(string) != 0) )
            h = (h + 1) & mask;
        
        return h;
    }
    
    void Insert(unsigned int where, unsigned int key, const String & string)
    {
        strings[where] = new String;
        *(strings[where]) = string;
        keys[where] = key;
        
        count++;
    }
};

class StringDoubleHash
{
protected:
    String        ** strings;
    double         * doubles;
    unsigned int      * keys;
    unsigned int count, size;
    unsigned int        mask;
    
public:
    StringDoubleHash(int startsize = 32);
    virtual ~StringDoubleHash();
    
    void Grow()    { SetSize(size * 2); }
    void Shrink()  { SetSize(size / 2); }
    
    void SetSize(int newsize);
    
    void Clear();
    
    int  Capacity() const { return size; }
    int  Entries() const  { return count; }
    
    double Double(int i) const { return doubles[i]; }
    double Double(const String & key) const
    {
        int index = Find(key);
        
        return index >= 0 ? doubles[index] : _NAN_;
    }
    
    void SetDouble(int i, double value)
    { doubles[i] = value; }
    void SetDouble(const String & key, double value)
    { Add(key, value); }
    
    int Add(const String & s, double value);
    int Find(const String & s, double defaultValue);
    int Find(const String & s) const;
    
    StringDoubleHash & operator = (const StringDoubleHash & rhs);
    
    const String & operator [] (int i) const { return *(strings[i]); }
    String & operator [] (int i) { return *(strings[i]); }
    //      String & String(int i) { return *(strings[i]); }
    
    void Delete(unsigned int index);
    void Delete(const String & key) { Delete(Find(key)); }
    
    bool SlotInUse(int index) { return strings[index] != NULL; }
    
private:
    
    unsigned int Iterate(unsigned int key, const String & string) const
    {
        unsigned int h = key & mask;
        
        while (  strings[h] != NULL &&
               ( keys[h] != key ||
                strings[h]->SlowCompare(string) != 0) )
            h = (h + 1) & mask;
        
        return h;
    }
    
    void Insert(unsigned int where, unsigned int key, const String & string)
    {
        strings[where] = new String;
        *(strings[where]) = string;
        keys[where] = key;
        
        count++;
    }
};

#endif
