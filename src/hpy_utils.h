#ifndef SRC_HPY_UTILS_H_
#define SRC_HPY_UTILS_H_

#include <hpy.h>

static inline int
HPy_SetAttrStringConstant_s(HPyContext *ctx, HPy obj, const char *name, const char *value)
{
    HPy h_value = HPyUnicode_FromString(ctx, value);
    if (HPy_IsNull(h_value)) {
        return -1;
    }
    int res = HPy_SetAttr_s(ctx, obj, name, h_value);
    HPy_Close(ctx, h_value);
    return res;
}

static inline int
HPy_SetAttrIntConstant_s(HPyContext *ctx, HPy obj, const char *name, long value)
{
    HPy h_value = HPyLong_FromLong(ctx, value);
    if (HPy_IsNull(h_value)) {
        return -1;
    }
    int res = HPy_SetAttr_s(ctx, obj, name, h_value);
    HPy_Close(ctx, h_value);
    return res;
}

static inline long
HPy_GetLongItem_i(HPyContext *ctx, HPy obj, HPy_ssize_t idx)
{
    long res;
    HPy h_item = HPy_GetItem_i(ctx, obj, idx);
    if (HPy_IsNull(h_item)) {
        return -1;
    }
    res = HPyLong_AsLong(ctx, h_item);
    HPy_Close(ctx, h_item);
    return res;
}

static inline int
HPyLong_Check(HPyContext *ctx, HPy obj)
{
    return HPy_TypeCheck(ctx, obj, ctx->h_LongType);
}

static inline int
HPyFloat_Check(HPyContext *ctx, HPy obj)
{
    return HPy_TypeCheck(ctx, obj, ctx->h_FloatType);
}

static inline HPy
HPy_CallMethod_s(HPyContext *ctx, const char *name, const HPy *args, size_t nargs, HPy kwnames)
{
    HPy h_name = HPyUnicode_FromString(ctx, name);
    if (HPy_IsNull(h_name)) {
        return HPy_NULL;
    }
    HPy h_result = HPy_CallMethod(ctx, h_name, args, nargs, kwnames);
    HPy_Close(ctx, h_name);
    return h_result;
}

#define HPyErr_BadInternalCall(ctx) HPyErr_Format(ctx, ctx->h_SystemError, "%s:%d: bad argument to internal function", __FILE__, __LINE__)

#endif /* SRC_HPY_UTILS_H_ */
