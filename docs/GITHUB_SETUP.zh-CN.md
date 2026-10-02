# 将准备好的代码展示在 GitHub

建议名称：`hare-simulation`

建议简介：`Simulation code for history-assisted residual estimation of long-term treatment effects with limited target outcomes.`

本文件是发布操作说明；准备工作没有创建远程仓库，也没有上传文件。

## 网页上传

1. 登录你希望使用的 GitHub 账号，打开 [创建仓库页面](https://github.com/new)。填写上述名称和简介。要对外展示时可选择 Public；Owner 选择你的账号。
2. 本地目录已经包含 README 和 `.gitignore`，创建空仓库即可。软件许可证需要由你决定，目前没有添加。
3. 解压交付的压缩包，将 `hare-simulation` **里面的文件和子目录**上传到仓库根目录，保留目录结构。不要只上传 ZIP 文件本身。
4. 使用 Upload files 或页面提示的上传入口，检查 README、`04_cloud_run.R`、`data/reference/`、`scripts/`、`figures/`、`tables/`、`docs/`、`.gitignore` 和 `.gitattributes` 均已包含。确认目录中不包含 `generated/` 本机环境报告。
5. 提交后检查仓库首页能显示 README 和结果图，并确认文件可浏览。此时仓库网址即可用于代码展示。

以上流程依据 GitHub 官方的 [创建仓库说明](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository) 和 [上传文件说明](https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository)。

## 推荐的展示表述

可以表述为：本仓库提供论文 simulation 的原始代码、保存结果，以及从保存结果重建图表和核对数值的辅助脚本。

暂时不应表述为：已在全新环境中逐位重现全部原始 Monte Carlo 结果。原云端依赖版本未保存，本次也没有重跑完整模拟。

## 将来要重跑时

按 README 安装依赖，先运行 `pilot 1`。它会运行一个预算、一次重复、三个 lambda 值。新输出保存在忽略目录 `04_lambda_transport_cloud/`；论文参考结果保存在 `data/reference/`。不要将不同代码版本的旧 checkpoint 混在一起。
