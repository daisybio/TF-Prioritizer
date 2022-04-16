import { NgModule } from '@angular/core';
import { BrowserModule } from '@angular/platform-browser';

import { AppComponent } from './app.component';
import { TfrankingEntryComponent } from './tfranking-entry/tfranking-entry.component';

@NgModule({
  declarations: [
    AppComponent,
    TfrankingEntryComponent
  ],
  imports: [
    BrowserModule
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
