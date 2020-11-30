#include <stdlib.h>
#include <string.h>

#include <stdio.h>

#include "smap.h"
#define SMAP_ARR_INITIAL_SIZE 4

typedef struct node {
    char* key;
    int* val;
    struct node* next;
} node;

struct _smap {
    int size;
    int arr_cap;
    int (*hash)(const char* s);
    node** arr;
};

int smap_default_hash(const char *s){
    //Horner's rule with x=31
    int h = s[0];
    for (int i = 1; i < strlen(s); i++) {
        h = 31*h + s[i];
    }
    if (h < 0) {
        h = -h;
    }
    if (h < 0) {
        h = 0;
    }
    //printf("PRODUCED HASH %i FOR \"%s\"\n", h, s);
    return h;
}

void smap_destroy_array(node** arr, int n, bool destroy_keys) {
    for (int i = 0; i < n; i++) {
        node* curr = arr[i];
        while (curr->next != NULL) {
            node* next = curr->next;
            if (destroy_keys) {
                free(curr->key);
            }
            free(curr);
            curr = next;
        }
        if (destroy_keys) {
            free(curr->key);
        }
        free(curr);
    }
    free(arr);
}

node* smap_create_empty_node() {
    node* new = malloc(sizeof(node));
    if (new == NULL) {
        return NULL;
    }
    new->key = NULL;
    new->val = NULL;
    new->next = NULL;
    return new;
}

node** smap_create_array(int cap) {
    node** arr = malloc(sizeof(node*) * cap);
    if (arr == NULL){
        return NULL;
    }
    for (int i = 0; i < cap; i++) {
        arr[i] = smap_create_empty_node();
        if (arr[i] == NULL) {
            smap_destroy_array(arr, SMAP_ARR_INITIAL_SIZE, 1);
            return NULL;
        }
    }
    return arr;
}

smap* smap_create(int (*h)(const char* s)){
    if (h == NULL) {
        return NULL;
    }
    smap* map = malloc(sizeof(smap));
    if (map == NULL) {
        return NULL;
    }
    map->arr = smap_create_array(SMAP_ARR_INITIAL_SIZE);
    if (map->arr == NULL){
        free(map);
        return NULL;
    }
    map->size = 0;
    map->arr_cap = SMAP_ARR_INITIAL_SIZE;
    map->hash = h;
    return map;
}

int smap_size(const smap *m) {
    return m != NULL? m->size : -1;
}

node* smap_find_key_node(const smap* m, const char* key) {
    int index = m->hash(key) % m->arr_cap;
    // printf("FINDING m[%s] AT m->arr[%i]\n", key, index);
    node* curr = m->arr[index];
    while (curr->next != NULL && strcmp(curr->key, key) != 0) {
        curr = curr->next;
    }
    if (curr->key == NULL) {
        // printf("   not found\n");
    }
    else {
    // printf("   found\n");
    }
    return curr;
}

bool smap_contains_key(const smap *m, const char *key) {
    if (key == NULL || m == NULL) {
        return 0;
    }
    node* curr = smap_find_key_node(m, key);
    return curr->key != NULL;
}

bool smap_reallocate(smap* m) {
    int new_cap = m->arr_cap * 2;
    node** new_arr = smap_create_array(new_cap);
    if (new_arr == NULL) {
        return 0;
    }
    for (int i = 0; i < m->arr_cap; i++) {
        node* old_curr = m->arr[i];
        while(old_curr->next != NULL){
            int index = m->hash(old_curr->key) % new_cap;
            // printf("(%s=>%i) moved to m->arr[%i]\n", old_curr->key, *(old_curr->val), index);
            node* new_curr = new_arr[index];
            while (new_curr->next != NULL) {
                new_curr = new_curr->next;
            }
            new_curr->key = old_curr->key;
            new_curr->val = old_curr->val;
            new_curr->next = smap_create_empty_node();
            if (new_curr->next == NULL) {
                smap_destroy_array(new_arr, new_cap, 0);
                return 0;
            }
            old_curr = old_curr->next;
        }
    }
    smap_destroy_array(m->arr, m->arr_cap, 0);
    m->arr = new_arr;
    m->arr_cap = new_cap;
    return 1;
}

bool smap_put(smap *m, const char *key, int *value){
    // printf("PUTTING (%s=>%i)\n  ", key, *value);
    if (key == NULL || m == NULL) {
        return 0;
    }
    node* curr = smap_find_key_node(m, key);
    if (curr->key == NULL) {
        if (m->size+1 == m->arr_cap) {
            //reallocate array and rehash elements
            smap_reallocate(m);
            curr = smap_find_key_node(m, key);
        }
        curr->next = smap_create_empty_node();
        if (curr->next == NULL) {
            return 0;
        }
        char* new_key = malloc(sizeof(char) * (strlen(key)+1));
        if (new_key == NULL) {
            free(curr->next);
            curr->next = NULL;
            return 0;
        }
        strcpy(new_key, key);
        curr->key = new_key;
        m->size++;
    }
    curr->val = value;
    return 1;
}

int *smap_get(smap *m, const char *key) {
    if (key == NULL || m == NULL) {
        return NULL;
    }
    node* curr = smap_find_key_node(m, key);
    if (curr->key == NULL) {
        return NULL;
    }
    else {
        return curr->val;
    }
}

void smap_print(smap* m) {
    for (int i = 0; i < m->arr_cap; i++) {
        printf("m->arr[%i] = ", i);
        node* curr = m->arr[i];
        while(curr->next != NULL){
            printf("(%s=>%i) ", curr->key, *(curr->val));
            curr = curr->next;
        }
        printf("\n");
    }
}

int *smap_remove(smap *m, const char *key) {
    if (key == NULL || m == NULL) {
        return NULL;
    }
    int index = m->hash(key) % m->arr_cap;
    // printf("FINDING m[%s] AT m->arr[%i]\n", key, index);
    node* curr = m->arr[index];
    node** ptr_to_curr = &(m->arr[index]);
    while (curr->next != NULL && strcmp(curr->key, key) != 0) {
        ptr_to_curr = &(curr->next);
        curr = curr->next;
    }
    if (curr->key != NULL) {
        *ptr_to_curr = curr->next;
        int* ret = curr->val;
        free(curr->key);
        free(curr);
        m->size--;
        return ret;
    }
    else {
        return NULL;
    }
}

void smap_destroy(smap *m) {
    if(m != NULL) {
        smap_destroy_array(m->arr, m->arr_cap, 1);
        free(m);
    }
}

void smap_for_each(smap *m, void (*f)(const char *, int *, void *), void *arg) {
    if (m == NULL || f == NULL) {
        return;
    }
    for (int i = 0; i < m->arr_cap; i++) {
        node* curr = m->arr[i];
        while(curr->next != NULL){
            f(curr->key, curr->val, arg);
            curr = curr->next;
        }
    }
}