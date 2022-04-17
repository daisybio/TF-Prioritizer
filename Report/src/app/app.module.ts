import {NgModule} from '@angular/core';
import {BrowserModule} from '@angular/platform-browser';

import {AppRoutingModule} from "./app-routing";
import {AppComponent} from './app.component';
import {HeaderComponent} from './accordion/header/header.component';
import {TfGroupAccordionComponent} from './accordion/tf-group-accordion/tf-group-accordion.component';
import {SingleTfAccordionComponent} from './accordion/single-tf-accordion/single-tf-accordion.component';
import {TfRankingComponent} from './tf-ranking/tf-ranking.component';
import {DataHeaderComponent} from './accordion/single-tf-accordion/data-header/data-header.component';
import {GeneIDComponent} from './accordion/single-tf-accordion/gene-id/gene-id.component';
import {DataContentComponent} from './accordion/single-tf-accordion/data-content/data-content.component';
import {ValidationComponent} from './validation/validation.component';
import { DistributionComponent } from './distribution/distribution.component';
import { RegressionComponent } from './regression/regression.component';
import { DocumentationComponent } from './documentation/documentation.component';
import { ParametersComponent } from './parameters/parameters.component';

@NgModule({
  declarations: [
    AppComponent,
    HeaderComponent,
    TfGroupAccordionComponent,
    SingleTfAccordionComponent,
    TfRankingComponent,
    DataHeaderComponent,
    GeneIDComponent,
    DataContentComponent,
    ValidationComponent,
    DistributionComponent,
    RegressionComponent,
    DocumentationComponent,
    ParametersComponent,
  ],
  imports: [
    BrowserModule,
    AppRoutingModule
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule {
}
