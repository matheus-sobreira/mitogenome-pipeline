# Script PowerShell: Constrói todas as imagens Docker do pipeline
# Execute uma vez antes de rodar o Nextflow.
#
# Uso:
#   cd mitogenome-pipeline
#   .\build_images.ps1

$ErrorActionPreference = "Stop"

$images = @(
    @{ tag = "mitogenome-pipeline/sra-tools:1.0";   path = "docker\sra-tools"   },
    @{ tag = "mitogenome-pipeline/fastqc:1.0";       path = "docker\fastqc"       },
    @{ tag = "mitogenome-pipeline/trim-galore:1.0";  path = "docker\trim-galore"  },
    @{ tag = "mitogenome-pipeline/novoplasty:1.0";   path = "docker\novoplasty"   }
)

Write-Host ""
Write-Host "=== Construindo imagens Docker do mitogenome-pipeline ===" -ForegroundColor Cyan
Write-Host ""

foreach ($img in $images) {
    Write-Host ">> Construindo: $($img.tag)" -ForegroundColor Yellow
    docker build -t $img.tag $img.path
    if ($LASTEXITCODE -ne 0) {
        Write-Host "ERRO ao construir $($img.tag)" -ForegroundColor Red
        exit 1
    }
    Write-Host "   OK: $($img.tag)" -ForegroundColor Green
    Write-Host ""
}

Write-Host "=== Todas as imagens construídas com sucesso! ===" -ForegroundColor Green
Write-Host ""
Write-Host "Imagens disponíveis:" -ForegroundColor Cyan
docker images | Select-String "mitogenome-pipeline"
